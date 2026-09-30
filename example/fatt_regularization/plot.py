#
# Plot the true and inverted models and the model errors
#
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans'],
                     'font.size': 11, 'axes.titlesize': 11, 'axes.labelsize': 11,
                     'xtick.labelsize': 10, 'ytick.labelsize': 10, 'legend.fontsize': 9})

nx = nz = 101
dx = 0.02
niter = 60


def load(path):
    """Read a LATTE model file as an (nz, nx) array."""
    return np.fromfile(path, np.float32).reshape(nx, nz).T.astype(np.float64)


vp = load('model/vp.bin')

# Model error: root-mean-square difference from the true model, from the
# initial model (iteration 0) to the last iteration
runs = [('test_noreg', 'No regularization', '#1f5fa8', '-'),
        ('test_tgpv', 'TGpV', '#c23b22', '-'),
        ('test_tv', 'TV', '#2a9d5c', '-'),
        ('test_tgpv_delayed', 'TGpV from iteration 11', '#c23b22', '--')]
error = {}
for run, *_ in runs:
    models = ['model/vp_init.bin'] + [f'{run}/iteration_{k}/model/updated_vp.bin' for k in range(1, niter + 1)]
    error[run] = [np.sqrt(np.mean((load(m) - vp)**2)) for m in models]

print(f'{"Run":<26}{"Final error (m/s)":>19}{"Minimum error (m/s)":>21}{"Iteration":>11}')
for run, label, *_ in runs:
    e = error[run]
    print(f'{label:<26}{e[-1]:>19.1f}{min(e):>21.1f}{np.argmin(e):>11d}')

fig, axes = plt.subplots(2, 3, figsize=(6.8, 4.9), constrained_layout=True)

# True and final inverted models
extent = [-0.5*dx, (nx - 0.5)*dx, (nz - 0.5)*dx, -0.5*dx]
panels = [(vp, '(a) True model')]
for (run, label, *_), tag in zip(runs, 'bcde'):
    panels.append((load(f'{run}/iteration_{niter}/model/updated_vp.bin'), f'({tag}) {label}'))
for ax, (m, title) in zip(axes.ravel()[:5], panels):
    im = ax.imshow(m, extent=extent, cmap='viridis', vmin=2800, vmax=3200, interpolation='none')
    ax.set_title(title, loc='left')
    ax.set_xticks([0, 1, 2])
    ax.set_yticks([0, 1, 2])
for ax in axes[:, 0]:
    ax.set_ylabel('Depth (km)')
for ax in axes[1, :2]:
    ax.set_xlabel('Distance (km)')
for ax in (axes[0, 1], axes[0, 2], axes[1, 1]):
    ax.tick_params(labelleft=False)
cb = fig.colorbar(im, ax=axes[:, :2].ravel().tolist() + [axes[0, 2]],
                  shrink=0.6, aspect=25, pad=0.02, location='bottom')
cb.set_label('P-wave velocity (m/s)')

# Model errors
ax = axes[1, 2]
for run, label, color, style in runs:
    ax.plot(range(0, niter + 1), error[run], color=color, ls=style, lw=1.6, label=label)
ax.set_xlabel('Iteration')
ax.set_ylabel('Model error (m/s)')
ax.set_title('(f) Model error', loc='left')
ax.set_xlim(0, niter)
ax.set_ylim(65, 155)
ax.legend(frameon=False, loc='upper right')

Path('figures').mkdir(exist_ok=True)
fig.savefig('figures/checkerboard_regularization.pdf')
fig.savefig('figures/checkerboard_regularization.png', dpi=200)
