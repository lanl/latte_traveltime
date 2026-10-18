"""2D test case for adjoint_receiver = nearest and linear.

Grid 121 x 61, 10 m spacing. Background v = 1800 + 1.2 z with a slow Gaussian
anomaly; the "true" model adds two more anomalies for realistic residuals.

Shot 1, source at (100, 0), receivers:
  sweep_x: z = 300, x = 900 ... 920 m in 1 m steps (moves across two cells)
  sweep_z: x = 1000, z = 200 ... 220 m in 1 m steps
  pair_same: two receivers in one cell that round to the same grid point
  pair_adj: two receivers in one cell that round to neighboring grid points
  surface: 44 surface receivers every 23.7 m
Shot 2, source at (640, 0): a borehole at x = 703 m, z = 15 ... 546 m every 9 m.
"""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lattetools import write_array

nx, nz, d = 121, 61, 10.0
x = np.arange(nx)*d
z = np.arange(nz)*d
Z, X = np.meshgrid(z, x, indexing="ij")


def gauss(x0, z0, s):
    return np.exp(-((X - x0)**2 + (Z - z0)**2)/(2*s**2))


v = 1800 + 1.2*Z - 250*gauss(600, 250, 100)
vtrue = v - 100*gauss(700, 60, 60) - 150*gauss(680, 350, 50)
write_array(v, "model/v.bin")
write_array(vtrue, "model/vtrue.bin")


def write_shot(path, src, rec):
    lines = ["1", "1", f"{src[0]:.9g} 0 {src[1]:.9g} 0", str(len(rec))]
    lines += [f"{r[0]:.9g} 0 {r[1]:.9g} 1" for r in rec]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("\n".join(lines) + "\n")


s = np.arange(21, dtype=float)
groups = {
    "sweep_x": [(900 + a, 300.0) for a in s],
    "sweep_z": [(1000.0, 200 + a) for a in s],
    "pair_same": [(801.5, 402.0), (804.0, 404.5)],
    "pair_adj": [(804.0, 402.0), (806.0, 404.0)],
    "surface": [(160 + 23.7*k, 0.0) for k in range(44)],
}
shot1 = []
index = {}
for name, rec in groups.items():
    index[name] = list(range(len(shot1), len(shot1) + len(rec)))
    shot1 += rec
write_shot("geometry/shot1_all.txt", (100.0, 0.0), shot1)
for name, rec in groups.items():
    write_shot(f"geometry/shot1_{name}.txt", (100.0, 0.0), rec)

vsp = [(703.0, 15.0 + 9*k) for k in range(60)]
write_shot("geometry/shot2_all.txt", (640.0, 0.0), vsp)
write_shot("geometry/shot2_vsp.txt", (640.0, 0.0), vsp)
index["vsp"] = list(range(60))

meta = dict(nx=nx, nz=nz, d=d, groups={k: [list(p) for p in r] for k, r in groups.items()},
            vsp=[list(p) for p in vsp], index=index, src1=[100.0, 0.0], src2=[640.0, 0.0])
Path("meta.json").write_text(json.dumps(meta, indent=1))
print("shot 1:", len(shot1), "receivers; shot 2:", len(vsp), "receivers")
