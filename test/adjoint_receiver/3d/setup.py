"""3D test case for adjoint_receiver = nearest and linear.

Grid 41 x 41 x 31 (x, y, z), 10 m spacing, v = 1800 + 1.2 z with a slow
Gaussian anomaly. Source at (50, 60, 0). Receivers:
  sweep: 11 positions along a diagonal, p = (300, 280, 200) + s*(1, 0.6, 0.8),
         s = 0, 2, ..., 20 m, crossing cell faces in x, y and z
  pair_same: two receivers in one cell that round to the same grid point
  pair_adj: two receivers in one cell that round to neighboring grid points
"""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lattetools import write_array

nx, ny, nz, d = 41, 41, 31, 10.0
x, y, z = np.arange(nx)*d, np.arange(ny)*d, np.arange(nz)*d
Z, Y, X = np.meshgrid(z, y, x, indexing="ij")
v = 1800 + 1.2*Z - 200*np.exp(-((X - 200)**2 + (Y - 200)**2 + (Z - 120)**2)/(2*80.0**2))
write_array(v, "model/v.bin")

s = np.arange(0, 21, 2.0)
groups = {
    "sweep": [(300 + a, 280 + 0.6*a, 200 + 0.8*a) for a in s],
    "pair_same": [(311.0, 262.0, 181.0), (313.5, 263.5, 183.5)],
    "pair_adj": [(313.5, 262.0, 181.0), (316.0, 263.0, 182.0)],
}
src = (50.0, 60.0, 0.0)


def write_shot(path, rec):
    lines = ["1", "1", f"{src[0]:.9g} {src[1]:.9g} {src[2]:.9g} 0", str(len(rec))]
    lines += [f"{r[0]:.9g} {r[1]:.9g} {r[2]:.9g} 1" for r in rec]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("\n".join(lines) + "\n")


allrec, index = [], {}
for name, rec in groups.items():
    index[name] = list(range(len(allrec), len(allrec) + len(rec)))
    allrec += rec
    write_shot(f"geometry/{name}.txt", rec)
write_shot("geometry/all.txt", allrec)
for name in ["pair_same", "pair_adj"]:
    np.array([1e-3, -1e-3], dtype=np.float32).tofile(f"geometry/res_{name}.bin")
# FD window, 1-based inclusive grid indices (x, y, z)
window = [24, 34, 21, 31, 13, 23]
json.dump(dict(nx=nx, ny=ny, nz=nz, d=d, src=src, groups={k: [list(p) for p in r] for k, r in groups.items()},
               index=index, window=window), open("meta.json", "w"), indent=1)
print(len(allrec), "receivers")
