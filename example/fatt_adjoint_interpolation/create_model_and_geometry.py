"""Create the model and the geometry of the adjoint-receiver example.

The model is 1200 m x 600 m on a 10 m grid, v = 1800 + 1.2 z m/s with a slow
Gaussian anomaly. One source at (100, 0) m. The receivers:

  sweep: one receiver moved from (1000, 200) m to (1000, 210) m in 1 m
         steps, i.e., from a grid point through a whole cell to the next one;
         each position has a residual of +1 ms on its own
  pair:  two receivers in one cell that round to the same grid point, with
         residuals of +1 ms and -1 ms

setup.json holds the source, the receivers, their residuals and the windows
in which compare_gradients.py computes finite-difference gradients.
"""

import json
from pathlib import Path

import numpy as np

nx, nz, d = 121, 61, 10.0
x = np.arange(nx)*d
z = np.arange(nz)*d
Z, X = np.meshgrid(z, x, indexing="ij")
v = 1800.0 + 1.2*Z - 250.0*np.exp(-((X - 600.0)**2 + (Z - 250.0)**2)/(2*100.0**2))

Path("model").mkdir(exist_ok=True)
# LATTE's layout: z runs fastest
v.astype(np.float32).T.tofile("model/v.bin")

source = [100.0, 0.0]
receivers = [{"x": 1000.0, "z": 200.0 + k, "residual": 1e-3, "case": f"sweep_{k}"} for k in range(11)]
receivers += [{"x": 801.5, "z": 402.0, "residual": 1e-3, "case": "pair"},
              {"x": 804.0, "z": 404.5, "residual": -1e-3, "case": "pair"}]

# All receivers in one shot, for the traveltimes and the finite differences
Path("geometry").mkdir(exist_ok=True)
lines = ["1", "1", f"{source[0]:g} 0 {source[1]:g} 0", str(len(receivers))]
lines += [f"{r['x']:g} 0 {r['z']:g} 1" for r in receivers]
Path("geometry/shot_1_geometry.txt").write_text("\n".join(lines) + "\n")
Path("geometry/geometry.txt").write_text("shot_1_geometry.txt\n")

# Windows (grid indices from 0, x first) of the finite-difference gradients:
# the end of the ray near the moving receiver, and the pair
windows = {"sweep": [94, 101, 18, 23], "pair": [76, 82, 37, 42]}
json.dump({"nx": nx, "nz": nz, "d": d, "source": source, "receivers": receivers, "windows": windows},
          open("setup.json", "w"), indent=1)
print(f"{len(receivers)} receivers; windows: " + ", ".join(
    f"{k} {(w[1] - w[0] + 1)*(w[3] - w[2] + 1)} grid points" for k, w in windows.items()))
