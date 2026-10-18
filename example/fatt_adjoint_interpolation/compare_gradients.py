"""Compute the gradients of the adjoint-receiver example with both receiver
treatments, and the finite-difference reference.

For each case in setup.json (each position of the moving receiver, and the
receiver pair), x_fatt2 computes the gradient of the misfit 1/2 sum_r e_r^2
at the model with adjoint_receiver = nearest and with adjoint_receiver =
linear. The observed traveltimes are the synthetic ones minus the residuals
of setup.json, so the residuals e_r are exactly those values. The reference
perturbs the velocity at each grid point of the windows by +-1% and runs
x_eikonal2, which gives the derivative of every receiver traveltime with
respect to the velocity there by central differences.

Usage: python3 compare_gradients.py [--bindir DIR] [--jobs N]
"""

import argparse
import json
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DELTA = 0.01


def run(bindir, program, param, cwd):
    """Run one LATTE executable on one rank; stop when it does not complete."""
    log = Path(cwd)/f"log_{Path(param).stem}.txt"
    with open(log, "w") as f:
        subprocess.run(["mpirun", "-np", "1", str(Path(bindir)/f"x_{program}2"), str(param)], cwd=cwd,
                       stdout=f, stderr=subprocess.STDOUT, env=dict(os.environ, OMP_NUM_THREADS="1"))
    if f"{program.upper()} Completed" not in log.read_text(errors="replace"):
        raise RuntimeError(f"x_{program}2 {param} in {cwd} did not complete; see {log}")


def write_shot(directory, source, receivers):
    directory.mkdir(parents=True, exist_ok=True)
    lines = ["1", "1", f"{source[0]:g} 0 {source[1]:g} 0", str(len(receivers))]
    lines += [f"{r['x']:g} 0 {r['z']:g} 1" for r in receivers]
    (directory/"shot_1_geometry.txt").write_text("\n".join(lines) + "\n")
    (directory/"geometry.txt").write_text("shot_1_geometry.txt\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--bindir", default=str(HERE/"../../bin"))
    parser.add_argument("--jobs", type=int, default=4, help="LATTE runs at the same time")
    args = parser.parse_args()
    bindir = Path(args.bindir).resolve()
    os.chdir(HERE)

    setup = json.load(open("setup.json"))
    nx, nz = setup["nx"], setup["nz"]
    receivers = setup["receivers"]
    v = np.fromfile("model/v.bin", dtype=np.float32).reshape(nx, nz).T
    eikonal = Path("param_eikonal.rb").read_text()
    gradient = Path("param_gradient.rb").read_text()
    work = Path("runs")
    if work.exists():
        shutil.rmtree(work)
    Path("results").mkdir(exist_ok=True)

    # The receiver traveltimes in the model
    base = work/"base"
    base.mkdir(parents=True)
    (base/"param.rb").write_text(eikonal + "file_geometry = ../../geometry/geometry.txt\nfile_vp = ../../model/v.bin\n")
    run(bindir, "eikonal", "param.rb", base)
    t = np.fromfile(base/"data/shot_1_traveltime_p.bin", dtype=np.float32)

    # The gradients of both settings, one x_fatt2 run per case and setting
    cases = list(dict.fromkeys(r["case"] for r in receivers))
    jobs = []
    for case in cases:
        index = [i for i, r in enumerate(receivers) if r["case"] == case]
        directory = work/case
        write_shot(directory/"geometry", setup["source"], [receivers[i] for i in index])
        (directory/"record").mkdir(parents=True)
        observed = t[index] - np.array([receivers[i]["residual"] for i in index], dtype=np.float32)
        observed.astype(np.float32).tofile(directory/"record/shot_1_traveltime_p.bin")
        for mode in ["nearest", "linear"]:
            (directory/f"param_{mode}.rb").write_text(
                gradient + f"file_geometry = geometry/geometry.txt\ndir_record = record\n"
                f"adjoint_receiver = {mode}\ndir_working = {mode}\n")
            jobs.append((case, mode, directory))
    with ThreadPoolExecutor(args.jobs) as pool:
        list(pool.map(lambda j: run(bindir, "fatt", f"param_{j[1]}.rb", j[2]), jobs))
    for case, mode, directory in jobs:
        g = np.fromfile(directory/mode/"iteration_1/model/grad_vp.bin", dtype=np.float32).reshape(nx, nz).T
        np.save(f"results/grad_{case}_{mode}.npy", g)

    # The finite-difference derivatives of all receiver traveltimes in the windows
    nodes = [(i, j) for w in setup["windows"].values() for i in range(w[0], w[1] + 1) for j in range(w[2], w[3] + 1)]

    def perturbed(node):
        i, j = node
        out = []
        for sign in (1, -1):
            directory = work/f"fd_{i}_{j}_{'p' if sign > 0 else 'm'}"
            (directory/"model").mkdir(parents=True)
            vv = v.copy()
            vv[j, i] *= 1 + sign*DELTA
            vv.astype(np.float32).T.tofile(directory/"model/v.bin")
            (directory/"param.rb").write_text(eikonal + "file_geometry = ../../geometry/geometry.txt\n")
            run(bindir, "eikonal", "param.rb", directory)
            out.append(np.fromfile(directory/"data/shot_1_traveltime_p.bin", dtype=np.float32).astype(np.float64))
            shutil.rmtree(directory)
        return (out[0] - out[1])/(2*DELTA*float(v[j, i]))

    with ThreadPoolExecutor(args.jobs) as pool:
        kernels = np.array(list(pool.map(perturbed, nodes)))
    np.save("results/fd_nodes.npy", np.array(nodes))
    np.save("results/fd_kernels.npy", kernels.T)
    np.save("results/traveltimes.npy", t)
    print(f"{len(jobs)} gradients and {2*len(nodes)} finite-difference runs in {work}")


if __name__ == "__main__":
    main()
