#
# Plot the true, initial and inverted models, a traveltime field, and the data misfit
#
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans'],
                     'font.size': 11, 'axes.titlesize': 11, 'axes.labelsize': 11,
                     'xtick.labelsize': 10, 'ytick.labelsize': 10})

nz, nx = 101, 301
d = 0.01
ns = 30
niter = 50
shot = 15


def load(path):
    """Read a LATTE model or traveltime field as an (nz, nx) array."""
    return np.fromfile(path, np.float32).reshape(nx, nz).T.astype(np.float64)


def traveltimes(directory):
    """Read the traveltimes of all shots."""
    return np.concatenate([np.fromfile(f'{directory}/shot_{i}_traveltime_p.bin', np.float32) for i in range(1, ns + 1)])


mask = load('model/mask.bin') == 1
vp = load('model/vp.bin')
models = [('True model', vp), ('Initial model', load('model/vp_init.bin')),
          (f'FATT, {niter} iterations', load(f'test/iteration_{niter}/model/updated_vp.bin'))]

# Model errors below the topography, and traveltime residuals
t = traveltimes('data')
for (label, m), it in zip(models[1:], (0, niter)):
    e = np.sqrt(np.mean((m[mask] - vp[mask])**2))
    r = np.sqrt(np.mean((traveltimes(f'test/iteration_{it}/synthetic') - t)**2))*1000
    print(f'{label:<24} RMS model error = {e:6.1f} m/s, RMS traveltime residual = {r:6.2f} ms')

# Source positions, from the fourth line of each shot file (x, y, z, t0)
src = np.array([open(f'geometry/shot_{i}_geometry.txt').read().split('\n')[3].split() for i in range(1, ns + 1)],
               dtype=float)

Path('figures').mkdir(exist_ok=True)

# Models; the region above the topography is blank
fig, axes = plt.subplots(3, 1, figsize=(6.8, 6.4), constrained_layout=True)
extent = [-0.5*d, (nx - 0.5)*d, (nz - 0.5)*d, -0.5*d]
for ax, (label, m), tag in zip(axes, models, 'abc'):
    im = ax.imshow(np.where(mask, m, np.nan), extent=extent, cmap='jet', vmin=2000, vmax=6000,
                   interpolation='none')
    ax.set_title(f'({tag}) {label}', loc='left')
    ax.set_ylabel('Depth (km)')
    ax.set_yticks([0, 0.5, 1])
for ax in axes[:2]:
    ax.tick_params(labelbottom=False)
axes[2].set_xlabel('Distance (km)')

# Traveltime field of one shot, every 0.05 s, and the sources
x, z = np.meshgrid(np.arange(nx)*d, np.arange(nz)*d)
tt = np.where(mask, load(f'snapshot/shot_{shot}_traveltime_p.bin'), np.nan)
axes[0].contour(x, z, tt, levels=np.arange(0.05, np.nanmax(tt), 0.05), colors='k', linewidths=0.8)
axes[0].plot(src[:, 0]/1000, src[:, 2]/1000, '*', color='w', mec='k', ms=9, clip_on=False)

cb = fig.colorbar(im, ax=axes, location='right', shrink=0.6, aspect=25)
cb.set_label('P-wave velocity (m/s)')
fig.savefig('figures/models.pdf')
fig.savefig('figures/models.png', dpi=200)

# Data misfit, normalized by its value in the initial model
misfit = np.loadtxt('test/data_misfit.txt')
fig, ax = plt.subplots(figsize=(3.5, 2.8), constrained_layout=True)
ax.semilogy(misfit[:, 0], misfit[:, 2], color='#1f5fa8', lw=1.8)
ax.set_xlabel('Iteration')
ax.set_ylabel('Normalized data misfit')
ax.set_xlim(0, niter)
ax.grid(True, which='both', color='0.85', lw=0.6)
fig.savefig('figures/misfit.pdf')
fig.savefig('figures/misfit.png', dpi=200)
