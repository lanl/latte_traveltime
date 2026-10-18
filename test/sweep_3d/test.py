"""Test the level ranges of the 3D fast sweeps on a grid with ny < nz.

A sweep of the 3D solvers visits the grid level by level. The range of j in
a level used nx + ny where it needs nx + nz, so on a grid with ny < nz every
sweep skipped some grid points. The other sweeps still reached them, but the
solver needed more iterations. With the default stopping settings
(sweep_stop_threshold = 1.0e-4, sweep_niter_max = 10), the traveltimes of
this grid differed from the fully converged ones by up to 256 us before the
fix, and by at most 0.5 us after it.

The test runs x_eikonal3 on a 41 x 17 x 31 grid at 20 m (ny < nz) with a
velocity gradient and a Gaussian anomaly, for 3 surface shots, once with the
default stopping settings and once converged fully. It passes when the two
traveltime fields agree within 10 us. When x_eikonal3_cuda exists, whose GPU
sweeps visit every point of a level, its fields must also equal those of
x_eikonal3 bit for bit.

python3 test.py [--bindir ../../bin] [--np 3]
"""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
nx, ny, nz, d = 41, 17, 31, 20.0
TOLERANCE = 10.0e-6


def setup():
    x, y, z = np.arange(nx)*d, np.arange(ny)*d, np.arange(nz)*d
    Z, Y, X = np.meshgrid(z, y, x, indexing="ij")
    v = 2000.0 + 1.5*Z - 300.0*np.exp(-((X - 420)**2 + (Y - 160)**2 + (Z - 300)**2)/(2*90.0**2))
    (HERE / "model").mkdir(exist_ok=True)
    # LATTE's layout: z runs fastest, then y, then x
    v.astype(np.float32).transpose(2, 1, 0).tofile(HERE / "model" / "vp.bin")
    receivers = [(rx, ry, 0.0) for rx in np.arange(0, 801, 40.0) for ry in np.arange(0, 321, 40.0)]
    (HERE / "geometry").mkdir(exist_ok=True)
    names = []
    for i, (sx, sy) in enumerate([(60.0, 60.0), (740.0, 160.0), (400.0, 260.0)], start=1):
        lines = [str(i), "1", f"{sx:g} {sy:g} 0 0", str(len(receivers))] + [f"{a:g} {b:g} {c:g} 1" for a, b, c in receivers]
        (HERE / "geometry" / f"shot_{i}_geometry.txt").write_text("\n".join(lines) + "\n")
        names.append(f"shot_{i}_geometry.txt")
    (HERE / "geometry" / "geometry.txt").write_text("\n".join(names) + "\n")


def run(executable, name, np_, converged=False):
    """x_eikonal on the model; the traveltime fields go to snap_<name>."""
    for folder in (f"data_{name}", f"snap_{name}"):
        shutil.rmtree(HERE / folder, ignore_errors=True)
    text = (f"nx = {nx}\nny = {ny}\nnz = {nz}\ndx = {d:g}\ndy = {d:g}\ndz = {d:g}\nns = 3\n"
            "file_geometry = geometry/geometry.txt\nmodel_name = vp\nfile_vp = model/vp.bin\n"
            f"dir_synthetic = data_{name}\nsnaps = 1\ndir_snapshot = snap_{name}\n")
    if converged:
        text += "sweep_stop_threshold = 1.0e-12\nsweep_niter_max = 200\n"
    (HERE / f"param_{name}.rb").write_text(text)
    with open(HERE / f"{name}.log", "w") as f:
        r = subprocess.run(["mpirun", "-np", str(np_), str(executable), f"param_{name}.rb"], cwd=HERE, stdout=f,
                           stderr=subprocess.STDOUT, env=dict(os.environ, OMP_NUM_THREADS="1"))
    if r.returncode != 0:
        sys.exit(f"{executable.name} failed; see {HERE / f'{name}.log'}")
    return [np.fromfile(HERE / f"snap_{name}" / f"shot_{i}_traveltime_p.bin", dtype=np.float32) for i in (1, 2, 3)]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bindir", default=str(HERE / "../../bin"))
    ap.add_argument("--np", type=int, default=3)
    args = ap.parse_args()
    bindir = Path(args.bindir).resolve()
    setup()

    default = run(bindir / "x_eikonal3", "default", args.np)
    converged = run(bindir / "x_eikonal3", "converged", args.np, converged=True)
    error = max(float(np.max(np.abs(a.astype(np.float64) - b))) for a, b in zip(default, converged))
    passed = error <= TOLERANCE
    print(f"x_eikonal3, default stopping against converged: largest difference {error*1e6:.2f} us "
          f"(tolerance {TOLERANCE*1e6:.0f} us): {'PASS' if passed else 'FAIL'}")

    if (bindir / "x_eikonal3_cuda").is_file():
        cuda = run(bindir / "x_eikonal3_cuda", "cuda", args.np)
        log = (HERE / "cuda.log").read_text()
        gpu = " on GPU " in log and "GPU fast sweeping failed" not in log
        same = all(np.array_equal(a, b) for a, b in zip(default, cuda))
        print(f"x_eikonal3_cuda against x_eikonal3: {'identical' if same else 'different'} fields"
              f"{'' if gpu else ' (the CUDA run found no GPU)'}: {'PASS' if same and gpu else 'FAIL'}")
        passed = passed and same and gpu
    print("PASS" if passed else "FAIL")
    sys.exit(0 if passed else 1)


if __name__ == "__main__":
    main()
