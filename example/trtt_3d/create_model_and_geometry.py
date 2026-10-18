"""Create the models and the geometry of the 3D TRTT example.

The grid is 61 x 41 x 41 points (x, y, z) at 10 m. The starting model is
v = 1500 + z m/s, and the true model adds a Gaussian anomaly of +200 m/s at
120 m depth. One reflector lies at about 260 m depth, curved along x and y.
24 surface sources, 100 m apart, record 31 x 21 surface receivers, 20 m apart.
"""

from pathlib import Path

import numpy as np

nx, ny, nz, d = 61, 41, 41, 10.0
x, y, z = np.arange(nx)*d, np.arange(ny)*d, np.arange(nz)*d
Z, Y, X = np.meshgrid(z, y, x, indexing="ij")
vp_init = 1500.0 + 1.0*Z
vp_true = vp_init + 200.0*np.exp(-((X - 300.0)**2/100.0**2 + (Y - 200.0)**2/100.0**2 + (Z - 120.0)**2/50.0**2)/2)
depth = 260.0 + 30.0*np.sin(2*np.pi*x/600.0)[None, :] + 15.0*np.cos(2*np.pi*y/400.0)[:, None]

# Reflector model: label 1 at the grid point nearest to the reflector in each
# column, plus the offset of the reflector from that point in grid cells
refl = np.zeros((nz, ny, nx))
iz = np.floor(depth/d + 0.5).astype(int)
jj, ii = np.meshgrid(np.arange(ny), np.arange(nx), indexing="ij")
refl[iz, jj, ii] = 1.0 + depth/d - iz

model = Path("model")
model.mkdir(exist_ok=True)
# LATTE's layout: z runs fastest, then y, then x
for name, w in (("vp_init", vp_init), ("vp_true", vp_true), ("refl", refl), ("refl_zero", 0*refl)):
    w.astype(np.float32).transpose(2, 1, 0).tofile(model/f"{name}.bin")
np.save(model/"depth.npy", depth)

geometry = Path("geometry")
geometry.mkdir(exist_ok=True)
receivers = [(rx, ry) for rx in np.arange(0.0, 601.0, 20.0) for ry in np.arange(0.0, 401.0, 20.0)]
names = []
for i, (sx, sy) in enumerate([(sx, sy) for sx in np.arange(50.0, 600.0, 100.0) for sy in np.arange(50.0, 400.0, 100.0)],
                             start=1):
    lines = [str(i), "1", f"{sx:g} {sy:g} 0 0", str(len(receivers))] + [f"{rx:g} {ry:g} 0 1" for rx, ry in receivers]
    (geometry/f"shot_{i}_geometry.txt").write_text("\n".join(lines) + "\n")
    names.append(f"shot_{i}_geometry.txt")
(geometry/"geometry.txt").write_text("\n".join(names) + "\n")
print(f"{len(names)} sources, {len(receivers)} receivers each")
