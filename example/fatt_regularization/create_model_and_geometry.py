#
# Create the true and initial models and the source-receiver geometry
#
import numpy as np
from pathlib import Path

nx = nz = 101
d = 20.0

# True model: a 5 x 5 checkerboard of 400 m cells with 3000 +/- 150 m/s;
# the last row and column of grid points belong to the last cells
x = np.arange(nx)*d
z = np.arange(nz)*d
X, Z = np.meshgrid(x, z)
ix = np.minimum((X//400).astype(int), 4)
iz = np.minimum((Z//400).astype(int), 4)
vp = 3000.0 + 150.0*np.where((ix + iz) % 2 == 0, 1.0, -1.0)

# Initial model: 3000 m/s everywhere
vp_init = np.full((nz, nx), 3000.0)

# LATTE reads raw float32 arrays with z as the fastest axis
Path('model').mkdir(exist_ok=True)
vp.T.astype(np.float32).tofile('model/vp.bin')
vp_init.T.astype(np.float32).tofile('model/vp_init.bin')

# Receivers every 40 m along the four sides, and 10 sources on each side
L = (nx - 1)*d
rec = []
for s in np.arange(0.0, L, 40.0):
    rec += [(s, 0.0), (L, s), (L - s, L), (0.0, L - s)]
rec = sorted(set(rec))
side = np.arange(100.0, L, 200.0)
src = [(p, 0.0) for p in side] + [(p, L) for p in side] + [(0.0, p) for p in side] + [(L, p) for p in side]

# Each shot file lists the shot ID, the source point (x, y, z, t0) and the receivers (x, y, z, weight)
Path('geometry').mkdir(exist_ok=True)
with open('geometry/geometry.txt', 'w') as f:
    for i, (sx, sz) in enumerate(src, start=1):
        f.write(f'shot_{i}_geometry.txt\n')
        lines = [str(i), '', '1', f'{sx:.3f} 0 {sz:.3f} 0', '', str(len(rec))]
        lines += [f'{rx:.3f} 0 {rz:.3f} 1' for rx, rz in rec]
        Path(f'geometry/shot_{i}_geometry.txt').write_text('\n'.join(lines) + '\n')

print(f'{len(src)} sources, {len(rec)} receivers per source')
