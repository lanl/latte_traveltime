"""Create the models and the geometry of the 2D TRTT example with faults.

RGM builds a layered model with four normal faults on a grid of 401 x 161
points at 10 m. The three reflectors are the RGT horizons at the largest
velocity jumps of the central column, so they step across the faults. The
starting model is the linear fit of the laterally averaged true velocity. 20
surface sources, 200 m apart, record 201 surface receivers, 20 m apart.
"""

from pathlib import Path

import numpy as np
import rgm

nx, nz, d = 401, 161, 10.0

# The faulted model and its relative geologic time (RGT); rgm.rgm2 of the
# Python interface is the Fortran type rgm2_curved of RGM
p = rgm.rgm2(n1=nz, n2=nx, nf=4, nl=12, seed=2026)
p.yn_rgt = True
p.yn_regular_fault = True
p.dip = [60.0, 75.0]
p.disp = [6.0, 12.0]
p.vmin = 1800.0
p.vmax = 3600.0
p.refl_smooth = 40.0
p.refl_height = [0.0, 15.0]
p.generate()
vp_true, rgt = p.vp.astype(np.float64), p.rgt.astype(np.float64)


def horizon(t):
    """The depth of the RGT value t in each column, by linear interpolation."""
    z = np.full(nx, np.nan)
    for j in range(nx):
        c = rgt[:, j]
        i = np.searchsorted(c, t, side="right") - 1
        if 0 <= i < nz - 1 and c[i + 1] > c[i]:
            z[j] = d*(i + (t - c[i])/(c[i + 1] - c[i]))
    return z


# The reflectors: the horizons at the largest velocity jump of the central
# column in each of three depth windows
j = nx//2
jump = np.diff(vp_true[:, j])
k = np.arange(nz - 1)
depth = []
for zmin, zmax in ((300.0, 650.0), (650.0, 1000.0), (1000.0, 1450.0)):
    w = (k*d >= zmin) & (k*d <= zmax)
    i = k[w][np.argmax(jump[w])]
    depth.append(horizon(0.5*(rgt[i, j] + rgt[i + 1, j])))
depth = np.array(depth)

# Reflector model: label l at the grid point nearest to reflector l in each
# column, plus the offset of the reflector from that point in grid cells
refl = np.zeros((nz, nx))
for label, z in enumerate(depth, start=1):
    ok = np.isfinite(z)
    iz = np.floor(z[ok]/d + 0.5).astype(int)
    refl[iz, np.nonzero(ok)[0]] = label + z[ok]/d - iz

# The starting model: the linear fit of the laterally averaged velocity
zz = np.arange(nz)*d
a, b = np.polyfit(zz, vp_true.mean(axis=1), 1)
vp_init = np.repeat((a*zz + b)[:, None], nx, axis=1)

model = Path("model")
model.mkdir(exist_ok=True)
# LATTE's layout: z runs fastest
for name, w in (("vp_init", vp_init), ("vp_true", vp_true), ("refl", refl), ("refl_zero", 0*refl)):
    w.astype(np.float32).T.tofile(model/f"{name}.bin")
np.save(model/"depth.npy", depth)

geometry = Path("geometry")
geometry.mkdir(exist_ok=True)
receivers = np.arange(0.0, 4001.0, 20.0)
names = []
for i, sx in enumerate(np.arange(100.0, 4000.0, 200.0), start=1):
    lines = [str(i), "1", f"{sx:g} 0 0 0", str(len(receivers))] + [f"{rx:g} 0 0 1" for rx in receivers]
    (geometry/f"shot_{i}_geometry.txt").write_text("\n".join(lines) + "\n")
    names.append(f"shot_{i}_geometry.txt")
(geometry/"geometry.txt").write_text("\n".join(names) + "\n")
print(f"{len(names)} sources, {len(receivers)} receivers each; reflectors at "
      + ", ".join(f"{np.nanmin(z):.0f}-{np.nanmax(z):.0f} m" for z in depth))
