"""Compare the CPU executables of the 3D solvers with their CUDA versions.

Each pair, x_eikonal3 and x_eikonal3_cuda, x_fatt3 and x_fatt3_cuda,
x_trtt3 and x_trtt3_cuda, and x_tloc3 and x_tloc3_cuda, runs the same
parameter file. The script times both runs and compares all their output
files byte by byte. The model has 101 x 81 x 41 grid points at 12 m, with a
velocity gradient, a Gaussian anomaly and one sinusoidal reflector. Four
surface sources record 13 x 11 surface receivers; tloc locates 8 sources at
depth from 7 x 6 surface receivers. fatt, trtt and tloc run 2 iterations.

python3 compare.py [--bindir ../../bin] [--np 4] [--programs eikonal fatt trtt tloc]

The exit status is 1 when a pair differs or a CUDA run finds no GPU.
"""

import argparse
import filecmp
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
nx, ny, nz, d = 101, 81, 41, 12.0


def write_geometry(folder, sources, receivers):
    """LATTE geometry files: one file per shot, listed in geometry.txt."""
    folder.mkdir(parents=True, exist_ok=True)
    names = []
    for i, (sx, sy, sz) in enumerate(sources, start=1):
        lines = [str(i), "1", f"{sx:g} {sy:g} {sz:g} 0", str(len(receivers))]
        lines += [f"{rx:g} {ry:g} {rz:g} 1" for rx, ry, rz in receivers]
        (folder / f"shot_{i}_geometry.txt").write_text("\n".join(lines) + "\n")
        names.append(f"shot_{i}_geometry.txt")
    (folder / "geometry.txt").write_text("\n".join(names) + "\n")


def setup():
    """The models, the geometries and the parameter files."""
    model = HERE / "model"
    model.mkdir(exist_ok=True)
    x, y, z = np.arange(nx)*d, np.arange(ny)*d, np.arange(nz)*d
    Z, Y, X = np.meshgrid(z, y, x, indexing="ij")
    vinit = 1500.0 + 1.0*Z
    vtrue = vinit + 200.0*np.exp(-((X - 600)**2 + (Y - 480)**2)/(2*150.0**2) - (Z - 200)**2/(2*70.0**2))
    # LATTE's layout: z runs fastest, then y, then x
    vinit.astype(np.float32).transpose(2, 1, 0).tofile(model / "vp_init.bin")
    vtrue.astype(np.float32).transpose(2, 1, 0).tofile(model / "vp_true.bin")
    # One reflector: label 1 plus the offset from the nearest grid point, in grid cells
    depth = 330.0 + 30.0*np.sin(2*np.pi*x/1200.0)[None, :] + 15.0*np.cos(2*np.pi*y/960.0)[:, None]
    refl = np.zeros((nz, ny, nx), dtype=np.float32)
    k = np.floor(depth/d + 0.5).astype(int)
    for j in range(ny):
        for i in range(nx):
            refl[k[j, i], j, i] = 1.0 + depth[j, i]/d - k[j, i]
    refl.transpose(2, 1, 0).tofile(model / "refl.bin")
    np.zeros((nx, ny, nz), dtype=np.float32).tofile(model / "refl_zero.bin")

    receivers = [(rx, ry, 0.0) for rx in np.linspace(0, 1200, 13) for ry in np.linspace(0, 960, 11)]
    surface = [(sx, sy, 0.0) for sx in (300.0, 900.0) for sy in (240.0, 720.0)]
    write_geometry(HERE / "geometry", surface, receivers)
    # tloc: 8 sources at depth with sparse surface receivers, and perturbed starting positions
    deep = np.array([(sx, sy, 300.0) for sx in (300.0, 500.0, 700.0, 900.0) for sy in (300.0, 660.0)])
    sparse = [(rx, ry, 0.0) for rx in np.linspace(0, 1200, 7) for ry in np.linspace(0, 960, 6)]
    write_geometry(HERE / "geometry_deep", [tuple(s) for s in deep], sparse)
    start = deep + np.random.default_rng(1).uniform(-30.0, 30.0, deep.shape)
    for k, name in enumerate(("sx", "sy", "sz")):
        start[:, k].astype(np.float32).tofile(model / f"{name}_init.bin")
    np.zeros(len(deep), dtype=np.float32).tofile(model / "st0_init.bin")

    grid = f"nx = {nx}\nny = {ny}\nnz = {nz}\ndx = {d:g}\ndy = {d:g}\ndz = {d:g}\n"
    surface_shots = "ns = 4\nfile_geometry = geometry/geometry.txt\n"
    inversion = ("niter_max = 2\nprocess_grad = smooth\ngrad_smoothx = 36\ngrad_smoothy = 36\ngrad_smoothz = 24\n"
                 "step_max_vp = 50\nmin_vp = 1000\nmax_vp = 4000\n")
    return {
        # The observed data, with the CPU executable
        "obs_fa": grid + surface_shots + "model_name = vp\nfile_vp = model/vp_true.bin\ndir_synthetic = obs_fa\n",
        "obs_refl": grid + surface_shots + "model_name = vp, refl\nfile_vp = model/vp_true.bin\nfile_refl = model/refl.bin\n"
                                           "dir_synthetic = obs_refl\n",
        "obs_deep": grid + "ns = 8\nfile_geometry = geometry_deep/geometry.txt\nmodel_name = vp\nfile_vp = model/vp_true.bin\n"
                           "dir_synthetic = obs_deep\n",
        # The comparisons; OUT is the output directory of each run
        "eikonal": grid + surface_shots + "model_name = vp, refl\nfile_vp = model/vp_true.bin\nfile_refl = model/refl.bin\n"
                                          "dir_synthetic = OUT\n",
        "fatt": grid + surface_shots + "model_update = vp\nfile_vp = model/vp_init.bin\ndir_record = obs_fa\n" + inversion
                + "dir_working = OUT\n",
        "trtt": grid + surface_shots + "model_update = vp, refl\nnrefl = 1\nfile_vp = model/vp_init.bin\n"
                                       "file_refl = model/refl_zero.bin\ndir_record = obs_refl\n" + inversion
                + "offset_max_refl_depth = 600\ndir_working = OUT\n",
        "tloc": grid + "ns = 8\nfile_geometry = geometry_deep/geometry.txt\nmodel_update = sx, sy, sz, st0\nmodel_aux = vp\n"
                       "file_vp = model/vp_true.bin\nfile_sx = model/sx_init.bin\nfile_sy = model/sy_init.bin\n"
                       "file_sz = model/sz_init.bin\nfile_st0 = model/st0_init.bin\ndir_record = obs_deep\nniter_max = 2\n"
                       "step_max_sx = 20\nstep_max_sy = 20\nstep_max_sz = 20\nstep_max_st0 = 0.01\nmisfit_type = ad\n"
                       "dir_working = OUT\n",
    }


def run(executable, text, out, np_, log):
    """Run one executable on a parameter file whose output directory is out."""
    if out.exists():
        shutil.rmtree(out)
    param = HERE / f"param_{out.name}.rb"
    param.write_text(text.replace("OUT", out.name))
    t = time.time()
    with open(log, "w") as f:
        r = subprocess.run(["mpirun", "-np", str(np_), str(executable), param.name], cwd=HERE, stdout=f,
                           stderr=subprocess.STDOUT, env=dict(os.environ, OMP_NUM_THREADS="1"))
    if r.returncode != 0:
        sys.exit(f"{executable.name} failed; see {log}")
    return time.time() - t


def compare(a, b):
    """All output files of two runs, byte by byte; logs and parameter dumps do not count."""
    def files(root):
        return sorted(p.relative_to(root) for p in root.rglob("*") if p.is_file() and "parameters." not in p.name
                      and p.suffix != ".log")
    fa, fb = files(a), files(b)
    differ = [str(f) for f in fa if f not in fb or not filecmp.cmp(a / f, b / f, shallow=False)]
    return len(fa), differ + [str(f) for f in fb if f not in fa]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bindir", default=str(HERE / "../../bin"))
    ap.add_argument("--np", type=int, default=4)
    ap.add_argument("--programs", nargs="+", default=["eikonal", "fatt", "trtt", "tloc"])
    args = ap.parse_args()
    bindir = Path(args.bindir).resolve()
    for p in ["eikonal"] + args.programs:
        for suffix in ("3", "3_cuda"):
            if not (bindir / f"x_{p}{suffix}").is_file():
                sys.exit(f"{bindir}/x_{p}{suffix} not found; build LATTE with use_cuda=on")

    params = setup()
    logs = HERE / "logs"
    logs.mkdir(exist_ok=True)
    for name in ("obs_fa", "obs_refl", "obs_deep"):
        run(bindir / "x_eikonal3", params[name], HERE / name, args.np, logs / f"{name}.log")

    rows, failed = [], False
    for p in args.programs:
        t_cpu = run(bindir / f"x_{p}3", params[p], HERE / f"{p}_cpu", args.np, logs / f"{p}_cpu.log")
        t_gpu = run(bindir / f"x_{p}3_cuda", params[p], HERE / f"{p}_cuda", args.np, logs / f"{p}_cuda.log")
        log = (logs / f"{p}_cuda.log").read_text()
        gpu = " on GPU " in log and "GPU fast sweeping failed" not in log
        n, differ = compare(HERE / f"{p}_cpu", HERE / f"{p}_cuda")
        failed |= bool(differ) or not gpu
        rows.append((p, t_cpu, t_gpu, n, differ, gpu))
        print(f"x_{p}3: {t_cpu:6.1f} s   x_{p}3_cuda: {t_gpu:6.1f} s   speedup {t_cpu/t_gpu:5.1f}   "
              f"{n} files, {'identical' if not differ else f'{len(differ)} differ: ' + ', '.join(differ[:3])}"
              f"{'' if gpu else '   (the CUDA run found no GPU)'}", flush=True)

    lines = [f"{'program':10s} {'CPU (s)':>8s} {'CUDA (s)':>9s} {'speedup':>8s} {'files':>6s}  result"]
    for p, t_cpu, t_gpu, n, differ, gpu in rows:
        result = "identical" if not differ else f"{len(differ)} files differ"
        if not gpu:
            result += ", no GPU"
        lines.append(f"{p:10s} {t_cpu:8.1f} {t_gpu:9.1f} {t_cpu/t_gpu:8.1f} {n:6d}  {result}")
    (HERE / "results.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
