"""Create the models and the geometry of the 2D TRTT example.

The grid is 301 x 151 points at 10 m. The starting model is v = 1500 + 0.8 z
m/s, and the true model adds a Gaussian anomaly of +300 m/s at 0.4 km depth.
One sinusoidal reflector lies between 0.8 and 1.0 km. 15 surface sources,
200 m apart, record 151 surface receivers, 20 m apart.
"""

from pathlib import Path

import numpy as np

nx, nz, d = 301, 151, 10.0
x, z = np.arange(nx)*d, np.arange(nz)*d
Z, X = np.meshgrid(z, x, indexing="ij")
vp_init = 1500.0 + 0.8*Z
vp_true = vp_init + 300.0*np.exp(-((X - 1500.0)**2/300.0**2 + (Z - 400.0)**2/100.0**2)/2)
depth = 900.0 + 100.0*np.sin(2*np.pi*x/3000.0)

# Reflector model: label 1 at the grid point nearest to the reflector in each
# column, plus the offset of the reflector from that point in grid cells
refl = np.zeros((nz, nx))
iz = np.floor(depth/d + 0.5).astype(int)
refl[iz, np.arange(nx)] = 1.0 + depth/d - iz

model = Path("model")
model.mkdir(exist_ok=True)
# LATTE's layout: z runs fastest
for name, w in (("vp_init", vp_init), ("vp_true", vp_true), ("refl", refl), ("refl_zero", 0*refl)):
    w.astype(np.float32).T.tofile(model/f"{name}.bin")
np.save(model/"depth.npy", depth)

geometry = Path("geometry")
geometry.mkdir(exist_ok=True)
receivers = np.arange(0.0, 3001.0, 20.0)
names = []
for i, sx in enumerate(np.arange(100.0, 3000.0, 200.0), start=1):
    lines = [str(i), "1", f"{sx:g} 0 0 0", str(len(receivers))] + [f"{rx:g} 0 0 1" for rx in receivers]
    (geometry/f"shot_{i}_geometry.txt").write_text("\n".join(lines) + "\n")
    names.append(f"shot_{i}_geometry.txt")
(geometry/"geometry.txt").write_text("\n".join(names) + "\n")
print(f"{len(names)} sources, {len(receivers)} receivers each")
