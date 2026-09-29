#!/usr/bin/env python
"""A source of several points in 2D: two point sources and a fault.

A shot's source may hold any number of points, each with its own start
time t0; the eikonal then gives, everywhere, the first arrival from any
of them. Four scenarios in a velocity gradient, run with x_eikonal2:

  two_points          two point sources firing together
  two_points_delayed  the same, the second fired 0.4 s later
  fault_bilateral     a dipping fault rupturing both ways from its middle
  fault_updip         the same fault rupturing up-dip from its deep end

The fault is a line of point sources 10 m apart; each point's t0 is its
distance along the fault from the hypocenter over the rupture velocity.

Needs the LATTE binaries in <repository>/bin, mpirun on the PATH, and
numpy and matplotlib. Each scenario runs in run_2d/<scenario>/, and the
figures go next to this script:

  traveltime_fields_2d.png/.pdf   the four traveltime fields
  surface_arrivals_2d.png/.pdf    first arrivals along the surface

Usage: python multi_source_2d.py
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

os.environ.setdefault("MPLBACKEND", "Agg")

HERE = Path(__file__).resolve().parent
LATTE = HERE.parents[1] / "bin" / "x_eikonal2"
MPIRUN = shutil.which("mpirun")
RUNS = HERE / "run_2d"

NX, NZ, D = 801, 401, 10.0            # 8 km x 4 km at 10 m
V0, GRADIENT = 1800.0, 0.7            # vp = V0 + GRADIENT * depth
RUPTURE_VELOCITY = 2500.0             # m/s
TWO = [(2500.0, 2000.0), (5500.0, 1200.0)]
FAULT = ((3000.0, 3000.0), (5000.0, 1200.0))      # deep end, shallow end
DELAY = 0.4                           # s, the second point source

# Arial, or a metric-compatible substitute; text stays editable in PDFs
STYLE = {"font.family": "sans-serif",
         "font.sans-serif": ["Arial", "Liberation Sans", "Arimo",
                             "Helvetica", "Nimbus Sans", "DejaVu Sans"],
         "pdf.fonttype": 42, "ps.fonttype": 42,
         "font.size": 11, "axes.labelsize": 11, "axes.titlesize": 11,
         "xtick.labelsize": 10, "ytick.labelsize": 10,
         "legend.fontsize": 10}


def fault_points(hypocenter: str):
    """Points 10 m apart along the fault, and each one's rupture time."""
    (x0, z0), (x1, z1) = FAULT
    length = np.hypot(x1 - x0, z1 - z0)
    s = np.arange(0.0, length + 1e-6, D)            # distance up-dip
    x = x0 + (x1 - x0) * s / length
    z = z0 + (z1 - z0) * s / length
    start = 0.5 * length if hypocenter == "middle" else 0.0
    t0 = np.abs(s - start) / RUPTURE_VELOCITY
    return np.column_stack([x, z, t0]), (x0 + (x1 - x0) * start / length,
                                         z0 + (z1 - z0) * start / length)


SCENARIOS = {
    "two_points": lambda: (np.array([[*TWO[0], 0.0], [*TWO[1], 0.0]]), None),
    "two_points_delayed": lambda: (np.array([[*TWO[0], 0.0],
                                             [*TWO[1], DELAY]]), None),
    "fault_bilateral": lambda: fault_points("middle"),
    "fault_updip": lambda: fault_points("deep end"),
}


def run(name: str) -> dict:
    """One shot whose source is all of the scenario's points."""
    folder = RUNS / name
    (folder / "model").mkdir(parents=True, exist_ok=True)
    (folder / "geometry").mkdir(exist_ok=True)
    depth = np.arange(NZ) * D
    vp = np.repeat((V0 + GRADIENT * depth)[:, None], NX, axis=1)
    np.ascontiguousarray(vp.T, dtype=np.float32).tofile(
        folder / "model" / "vp.bin")                    # z fastest
    points, hypocenter = SCENARIOS[name]()
    receivers = np.arange(0.0, (NX - 1) * D + 1e-6, 50.0)
    lines = ["1", "", str(len(points))]
    lines += [f"{x:.3f} 0 {z:.3f} {t:.6f}" for x, z, t in points]
    lines += ["", str(receivers.size)]
    lines += [f"{x:.3f} 0 0 1" for x in receivers]
    (folder / "geometry" / "shot_1_geometry.txt").write_text(
        "\n".join(lines) + "\n")
    (folder / "geometry" / "geometry.txt").write_text("shot_1_geometry.txt\n")
    (folder / "param_eikonal.rb").write_text("\n".join([
        f"nx = {NX}", f"nz = {NZ}", f"dx = {D:g}", f"dz = {D:g}", "ns = 1",
        "file_geometry = geometry/geometry.txt",
        "which_medium = acoustic-iso", "model_name = vp",
        "file_vp = model/vp.bin", "dir_synthetic = synthetic",
        "snaps = 1", "dir_snapshot = snapshot"]) + "\n")
    t_start = time.time()
    env = dict(os.environ)
    env.setdefault("OMP_NUM_THREADS", str(os.cpu_count() or 4))
    with open(folder / "eikonal.log", "w") as log:
        code = subprocess.run([MPIRUN, "-np", "1", str(LATTE),
                               "param_eikonal.rb"], cwd=folder, stdout=log,
                              stderr=subprocess.STDOUT, env=env).returncode
    if code != 0:
        raise RuntimeError(f"{name}: x_eikonal2 exited {code}; see "
                           f"{folder / 'eikonal.log'}")
    field = np.fromfile(folder / "snapshot" / "shot_1_traveltime_p.bin",
                        np.float32).reshape(NX, NZ).T
    arrivals = np.fromfile(folder / "synthetic" / "shot_1_traveltime_p.bin",
                           np.float32)
    return {"field": field, "arrivals": arrivals, "receivers": receivers,
            "points": points, "hypocenter": hypocenter,
            "seconds": time.time() - t_start}


def figures(results: dict) -> None:
    import matplotlib.pyplot as plt

    plt.rcParams.update(STYLE)
    km = 1e-3
    extent = [0.0, (NX - 1) * D * km, (NZ - 1) * D * km, 0.0]
    vmax = max(float(r["field"].max()) for r in results.values())
    levels = np.arange(0.1, vmax, 0.1)

    titles = {"two_points": "Two sources, fired together",
              "two_points_delayed": f"Right source fired {DELAY:g} s later",
              "fault_bilateral": "Fault, bilateral rupture",
              "fault_updip": "Fault, up-dip rupture"}
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 3.55), sharex=True,
                             sharey=True, constrained_layout=True)
    for ax, (name, r), letter in zip(axes.ravel(), results.items(), "abcd"):
        image = ax.imshow(r["field"], extent=extent, cmap="viridis",
                          vmin=0.0, vmax=vmax, aspect="equal",
                          interpolation="none")
        xs = (np.arange(NX) * D) * km
        zs = (np.arange(NZ) * D) * km
        ax.contour(xs, zs, r["field"], levels=levels, colors="white",
                   linewidths=0.6)
        pts = r["points"]
        if name.startswith("fault"):
            ax.plot(pts[:, 0] * km, pts[:, 1] * km, color="black", lw=2.2)
            hx, hz = r["hypocenter"]
            ax.plot(hx * km, hz * km, marker="*", color="#e8412c", ms=13,
                    mec="black", mew=0.8)
        else:
            ax.plot(pts[:, 0] * km, pts[:, 1] * km, "o", color="#e8412c",
                    ms=7, mec="black", mew=0.8)
        ax.set_title(f"({letter}) {titles[name]}", loc="left")
    for ax in axes[1, :]:
        ax.set_xlabel("Distance (km)")
    for ax in axes[:, 0]:
        ax.set_ylabel("Depth (km)")
    bar = fig.colorbar(image, ax=axes, shrink=0.95, pad=0.015, aspect=30)
    bar.set_label("Traveltime (s)")
    for suffix in ("png", "pdf"):
        fig.savefig(HERE / f"traveltime_fields_2d.{suffix}", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.0, 3.0), constrained_layout=True)
    styles = {"two_points": ("#1f5fa8", "-"),
              "two_points_delayed": ("#1f5fa8", "--"),
              "fault_bilateral": ("#c23b22", "-"),
              "fault_updip": ("#c23b22", "--")}
    for name, r in results.items():
        color, dash = styles[name]
        ax.plot(r["receivers"] * km, r["arrivals"], color=color,
                linestyle=dash, lw=1.8, label=titles[name])
    ax.set_xlabel("Distance (km)")
    ax.set_ylabel("First arrival (s)")
    ax.set_xlim(0.0, (NX - 1) * D * km)
    ax.invert_yaxis()
    ax.legend(loc="lower center", frameon=False)
    for suffix in ("png", "pdf"):
        fig.savefig(HERE / f"surface_arrivals_2d.{suffix}", dpi=300)
    plt.close(fig)


def main() -> int:
    if MPIRUN is None or not LATTE.exists():
        sys.exit(f"needs mpirun on the PATH and {LATTE}; build LATTE first")
    results = {}
    for name in SCENARIOS:
        results[name] = run(name)
        print(f"{name}: {len(results[name]['points'])} source point(s), "
              f"{results[name]['seconds']:.1f} s", flush=True)
    figures(results)
    print(f"figures in {HERE}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
