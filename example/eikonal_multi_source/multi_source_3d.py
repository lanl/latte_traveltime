#!/usr/bin/env python
"""A source of several points in 3D: two point sources and a rupturing
fault plane.

The 3D companion of multi_source_2d.py. A 6 x 6 x 3 km velocity
gradient at 25 m, run with x_eikonal3:

  two_points  two point sources on the plane y = 3 km, firing together
  fault       a rectangular fault, 3 km along strike (y) and 2 km down a
              45-degree dip to the east, as a grid of point sources 25 m
              apart; the rupture spreads over the plane from the center
              of its bottom edge at 2.5 km/s, which sets each point's t0

Needs the LATTE binaries in <repository>/bin, mpirun on the PATH, and
numpy and matplotlib. Each scenario runs in run_3d/<scenario>/, and the
figure traveltime_fields_3d.png/.pdf goes next to this script: for each
scenario, the first arrivals on the surface and on the vertical section
y = 3 km -- through both point sources, and across the fault's strike
through its hypocenter.

Usage: python multi_source_3d.py
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
LATTE = HERE.parents[1] / "bin" / "x_eikonal3"
MPIRUN = shutil.which("mpirun")
RUNS = HERE / "run_3d"

NX, NY, NZ, D = 241, 241, 121, 25.0   # 6 x 6 x 3 km at 25 m
V0, GRADIENT = 1800.0, 0.7            # vp = V0 + GRADIENT * depth
SECTION_Y = 3000.0                    # the vertical section shown
TWO = [(2000.0, SECTION_Y, 1500.0), (4500.0, SECTION_Y, 900.0)]

# the fault: strike along y, dipping 45 degrees to the east (+x)
STRIKE = (1500.0, 4500.0)             # y of its two ends
TOP_X, TOP_Z = 2500.0, 800.0          # its top edge
WIDTH, DIP = 2000.0, np.radians(45.0)
RUPTURE_VELOCITY = 2500.0             # m/s

# Arial, or a metric-compatible substitute; text stays editable in PDFs
STYLE = {"font.family": "sans-serif",
         "font.sans-serif": ["Arial", "Liberation Sans", "Arimo",
                             "Helvetica", "Nimbus Sans", "DejaVu Sans"],
         "pdf.fonttype": 42, "ps.fonttype": 42,
         "font.size": 11, "axes.labelsize": 11, "axes.titlesize": 11,
         "xtick.labelsize": 10, "ytick.labelsize": 10}


def fault_points():
    """A grid of point sources over the plane, and each one's rupture
    time from the center of the bottom edge."""
    s = np.arange(0.0, STRIKE[1] - STRIKE[0] + 1e-6, D)      # along strike
    w = np.arange(0.0, WIDTH + 1e-6, D)                       # down dip
    ss, ww = np.meshgrid(s, w, indexing="ij")
    x = TOP_X + ww * np.cos(DIP)
    y = STRIKE[0] + ss
    z = TOP_Z + ww * np.sin(DIP)
    s_h, w_h = 0.5 * (STRIKE[1] - STRIKE[0]), WIDTH
    t0 = np.hypot(ss - s_h, ww - w_h) / RUPTURE_VELOCITY
    hypocenter = (TOP_X + w_h * np.cos(DIP), STRIKE[0] + s_h,
                  TOP_Z + w_h * np.sin(DIP))
    return np.column_stack([x.ravel(), y.ravel(), z.ravel(),
                            t0.ravel()]), hypocenter


SCENARIOS = {
    "two_points": lambda: (np.array([[*TWO[0], 0.0], [*TWO[1], 0.0]]), None),
    "fault": fault_points,
}


def run(name: str) -> dict:
    folder = RUNS / name
    (folder / "model").mkdir(parents=True, exist_ok=True)
    (folder / "geometry").mkdir(exist_ok=True)
    depth = np.arange(NZ) * D
    vp = np.broadcast_to((V0 + GRADIENT * depth)[:, None, None],
                         (NZ, NY, NX))
    np.ascontiguousarray(vp.T, dtype=np.float32).tofile(
        folder / "model" / "vp.bin")                    # z fastest
    points, hypocenter = SCENARIOS[name]()
    grid = np.arange(0.0, (NX - 1) * D + 1e-6, 250.0)
    rx, ry = np.meshgrid(grid, grid, indexing="ij")
    lines = ["1", "", str(len(points))]
    lines += [f"{x:.3f} {y:.3f} {z:.3f} {t:.6f}" for x, y, z, t in points]
    lines += ["", str(rx.size)]
    lines += [f"{x:.3f} {y:.3f} 0 1" for x, y in zip(rx.ravel(), ry.ravel())]
    (folder / "geometry" / "shot_1_geometry.txt").write_text(
        "\n".join(lines) + "\n")
    (folder / "geometry" / "geometry.txt").write_text("shot_1_geometry.txt\n")
    (folder / "param_eikonal.rb").write_text("\n".join([
        f"nx = {NX}", f"ny = {NY}", f"nz = {NZ}", f"dx = {D:g}",
        f"dy = {D:g}", f"dz = {D:g}", "ns = 1",
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
        raise RuntimeError(f"{name}: x_eikonal3 exited {code}; see "
                           f"{folder / 'eikonal.log'}")
    field = np.fromfile(folder / "snapshot" / "shot_1_traveltime_p.bin",
                        np.float32).reshape(NX, NY, NZ).transpose(2, 1, 0)
    return {"field": field, "points": points, "hypocenter": hypocenter,
            "seconds": time.time() - t_start}


def figure(results: dict) -> None:
    import matplotlib.pyplot as plt

    plt.rcParams.update(STYLE)
    km = 1e-3
    xs = np.arange(NX) * D * km
    ys = np.arange(NY) * D * km
    zs = np.arange(NZ) * D * km
    j = int(round(SECTION_Y / D))
    vmax = max(float(max(r["field"][0].max(), r["field"][:, j, :].max()))
               for r in results.values())
    levels = np.arange(0.1, vmax, 0.1)
    style = dict(cmap="viridis", vmin=0.0, vmax=vmax, interpolation="none")

    fig = plt.figure(figsize=(7.0, 5.0), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, width_ratios=[1, 2])
    panels = [[fig.add_subplot(grid[r, c]) for c in range(2)]
              for r in range(2)]
    titles = {"two_points": "Two sources", "fault": "Fault"}
    letters = iter("abcd")
    for row, (name, r) in zip(panels, results.items()):
        surface, section = r["field"][0], r["field"][:, j, :]
        pts = r["points"]

        ax = row[0]
        image = ax.imshow(surface, extent=[xs[0], xs[-1], ys[0], ys[-1]],
                          origin="lower", aspect="equal", **style)
        ax.contour(xs, ys, surface, levels=levels, colors="white",
                   linewidths=0.6)
        if name == "fault":
            corners = np.array([[TOP_X, STRIKE[0]], [TOP_X, STRIKE[1]],
                                [pts[:, 0].max(), STRIKE[1]],
                                [pts[:, 0].max(), STRIKE[0]],
                                [TOP_X, STRIKE[0]]]) * km
            ax.plot(corners[:, 0], corners[:, 1], color="black", lw=1.4,
                    ls="--")
            hx, hy, _hz = r["hypocenter"]
            ax.plot(hx * km, hy * km, marker="*", color="#e8412c", ms=13,
                    mec="black", mew=0.8)
        else:
            ax.plot(pts[:, 0] * km, pts[:, 1] * km, "o", color="#e8412c",
                    ms=7, mec="black", mew=0.8)
        ax.axhline(SECTION_Y * km, color="white", lw=1.0, ls=":")
        ax.set_title(f"({next(letters)}) {titles[name]}, surface", loc="left")
        ax.set_xlabel("x (km)")
        ax.set_ylabel("y (km)")

        ax = row[1]
        ax.imshow(section, extent=[xs[0], xs[-1], zs[-1], zs[0]],
                  aspect="equal", **style)
        ax.contour(xs, zs, section, levels=levels, colors="white",
                   linewidths=0.6)
        if name == "fault":
            on = np.abs(pts[:, 1] - SECTION_Y) < 0.5 * D
            ax.plot(pts[on, 0] * km, pts[on, 2] * km, color="black", lw=2.2)
            hx, _hy, hz = r["hypocenter"]
            ax.plot(hx * km, hz * km, marker="*", color="#e8412c", ms=13,
                    mec="black", mew=0.8)
        else:
            ax.plot(pts[:, 0] * km, pts[:, 2] * km, "o", color="#e8412c",
                    ms=7, mec="black", mew=0.8)
        ax.set_title(f"({next(letters)}) {titles[name]}, section "
                     f"y = {SECTION_Y * km:g} km", loc="left")
        ax.set_xlabel("x (km)")
        ax.set_ylabel("Depth (km)")

    bar = fig.colorbar(image, ax=[a for row in panels for a in row],
                       shrink=0.9, pad=0.015, aspect=30)
    bar.set_label("Traveltime (s)")
    for suffix in ("png", "pdf"):
        fig.savefig(HERE / f"traveltime_fields_3d.{suffix}", dpi=300)
    plt.close(fig)


def main() -> int:
    if MPIRUN is None or not LATTE.exists():
        sys.exit(f"needs mpirun on the PATH and {LATTE}; build LATTE first")
    results = {}
    for name in SCENARIOS:
        results[name] = run(name)
        print(f"{name}: {len(results[name]['points'])} source point(s), "
              f"{results[name]['seconds']:.1f} s", flush=True)
    figure(results)
    print(f"figure in {HERE}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
