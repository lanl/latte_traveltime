#
# Plot the true, initial and inverted models, and the data misfit
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

nz, ny, nx = 61, 201, 201
d = 0.01
ns, nr = 16, 1681
niter = 9


def load(path):
    """Read a LATTE model file as an (nz, ny, nx) array."""
    return np.fromfile(path, np.float32).reshape(nx, ny, nz).transpose(2, 1, 0).astype(np.float64)


def traveltimes(directory):
    """Read the traveltimes of all shots."""
    return np.concatenate([np.fromfile(f'{directory}/shot_{i}_traveltime_p.bin', np.float32) for i in range(1, ns + 1)])


vp = load('model/vp.bin')
models = [('True', vp), ('Initial', load('model/vp_init.bin')),
          ('FATT', load(f'test/iteration_{niter}/model/updated_vp.bin'))]

# Model errors and traveltime residuals of the initial and inverted models
t = traveltimes('data')
for (label, m), it in zip(models[1:], (0, niter)):
    e = np.sqrt(np.mean((m - vp)**2))
    r = np.sqrt(np.mean((traveltimes(f'test/iteration_{it}/synthetic') - t)**2))*1000
    print(f'{label:<8} model: RMS model error = {e:6.1f} m/s, RMS traveltime residual = {r:6.2f} ms')

Path('figures').mkdir(exist_ok=True)

# Vertical sections through the center of the model
fig, axes = plt.subplots(3, 2, figsize=(7.0, 4.6), constrained_layout=True)
extent = [-0.5*d, (nx - 0.5)*d, (nz - 0.5)*d, -0.5*d]
tag = iter('abcdef')
for row, (label, m) in zip(axes, models):
    for ax, w, axis in zip(row, (m[:, ny//2, :], m[:, :, nx//2]), ('y', 'x')):
        im = ax.imshow(w, extent=extent, cmap='jet', vmin=1000, vmax=3000, interpolation='none')
        ax.set_title(f'({next(tag)}) {label}, {axis} = 1 km', loc='left')
        ax.set_yticks([0, 0.3, 0.6])
for ax in axes[:, 0]:
    ax.set_ylabel('Depth (km)')
for ax in axes[:, 1]:
    ax.tick_params(labelleft=False)
for ax in axes[:2, :].ravel():
    ax.tick_params(labelbottom=False)
axes[2, 0].set_xlabel('X (km)')
axes[2, 1].set_xlabel('Y (km)')
cb = fig.colorbar(im, ax=axes, location='bottom', shrink=0.5, aspect=30)
cb.set_label('P-wave velocity (m/s)')
fig.savefig('figures/sections.pdf')
fig.savefig('figures/sections.png', dpi=200)

# Depth slices of the true and inverted models; each depth has its own color range
depths = [0.1, 0.25, 0.4]
fig, axes = plt.subplots(2, 3, figsize=(7.0, 5.4), constrained_layout=True)
extent = [-0.5*d, (nx - 0.5)*d, -0.5*d, (ny - 0.5)*d]
tag = iter('abcdef')
for row, (label, m) in zip(axes, (models[0], models[2])):
    for ax, z in zip(row, depths):
        w = vp[round(z/d)]
        im = ax.imshow(m[round(z/d)], extent=extent, origin='lower', cmap='jet',
                       vmin=w.min(), vmax=w.max(), interpolation='none')
        ax.set_title(f'({next(tag)}) {label}, z = {z:g} km', loc='left')
        ax.set_xticks([0, 1, 2])
        ax.set_yticks([0, 1, 2])
for ax in axes[:, 0]:
    ax.set_ylabel('Y (km)')
for ax in axes[:, 1:].ravel():
    ax.tick_params(labelleft=False)
for ax in axes[0, :]:
    ax.tick_params(labelbottom=False)
for ax in axes[1, :]:
    ax.set_xlabel('X (km)')
for j in range(3):
    cb = fig.colorbar(axes[1, j].images[0], ax=axes[:, j], location='bottom', aspect=15)
    cb.set_label('Velocity (m/s)')
    cb.locator = matplotlib.ticker.MaxNLocator(4)
    cb.update_ticks()
fig.savefig('figures/depth_slices.pdf')
fig.savefig('figures/depth_slices.png', dpi=200)

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
