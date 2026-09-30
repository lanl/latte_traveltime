#
# Check the TGpV regularization against the FATT without regularization, and compare the models
#
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans'],
                     'font.size': 12, 'axes.titlesize': 12, 'axes.labelsize': 12,
                     'xtick.labelsize': 11, 'ytick.labelsize': 11, 'legend.fontsize': 11})

nz, nx = 101, 301
d = 0.01
ns = 30
niter = 50

# Output directory, label, regularization strength in each iteration, line color and style
runs = [('test', 'No regularization', lambda k: 0.0, '#222222', '-'),
        ('test_tgpv_0.1', 'TGpV, 0.1', lambda k: 0.1, '#1f77b4', '-'),
        ('test_tgpv_0.2', 'TGpV, 0.2', lambda k: 0.2, '#2ca02c', '-'),
        ('test_tgpv_0.5', 'TGpV, 0.5', lambda k: 0.5, '#d62728', '-'),
        ('test_tgpv_delayed', 'TGpV, 0.2 from iteration 11', lambda k: 0.2 if k >= 11 else 0.0, '#2ca02c', '--')]


def load(path):
    """Read a LATTE model or gradient as an (nz, nx) array."""
    return np.fromfile(path, np.float32).reshape(nx, nz).T.astype(np.float64)


def rms(a):
    return np.sqrt(np.mean(a**2))


def traveltimes(directory):
    """Read the traveltimes of all shots."""
    return np.concatenate([np.fromfile(f'{directory}/shot_{i}_traveltime_p.bin', np.float32) for i in range(1, ns + 1)])


def model(run, k):
    """Updated model of iteration k, which iteration k + 1 uses; iteration 0 gives the initial model."""
    return load(f'{run}/iteration_{k}/model/updated_vp.bin') if k > 0 else vp0


mask = load('model/mask.bin') == 1
vp, vp0 = load('model/vp.bin'), load('model/vp_init.bin')

# Checks of the regularization, against the run without regularization
print('Checks of the regularization')
for run, label, strength, *_ in runs[1:]:
    # LATTE denoises the updated model when the next iteration has a positive strength
    denoised = [k for k in range(1, niter + 1) if Path(f'{run}/iteration_{k}/model/reg_vp.bin').exists()]
    expected = [k for k in range(1, niter + 1) if strength(k + 1) > 0]
    # The first regularization term comes in the iteration after the first denoising,
    # and the earlier iterations equal those without regularization
    k0 = expected[0] + 1
    same = all(np.array_equal(np.fromfile(f'{run}/iteration_{k}/model/{f}_vp.bin', np.float32),
                              np.fromfile(f'test/iteration_{k}/model/{f}_vp.bin', np.float32))
               for k in range(1, k0) for f in ('grad', 'srch', 'updated'))
    # The term is coef*(m - mu) below the topography, where mu is the denoised model of the
    # previous iteration and coef = strength*RMS(gradient)/RMS(m - mu)
    g = load(f'test/iteration_{k0}/model/grad_vp.bin')
    r = model(run, k0 - 1) - load(f'{run}/iteration_{k0 - 1}/model/reg_vp.bin')
    term = load(f'{run}/iteration_{k0}/model/grad_vp.bin') - g
    expected_term = strength(k0)*rms(g)/rms(r)*r*mask
    mismatch = np.abs(term - expected_term).max()/np.abs(expected_term).max()
    # The velocity above the topography stays zero
    air = all(np.all(model(run, k)[~mask] == 0) for k in range(1, niter + 1))
    print(f'  {label}')
    print(f'    denoising in iterations {denoised[0]} to {denoised[-1]}: {"OK" if denoised == expected else "WRONG"}')
    print(f'    iterations 1 to {k0 - 1} equal those without regularization: {"OK" if same else "WRONG"}')
    print(f'    regularization term in iteration {k0}, relative mismatch {mismatch:.1e}: '
          f'{"OK" if mismatch < 1e-5 else "WRONG"}')
    print(f'    zero velocity above the topography: {"OK" if air else "WRONG"}')

# Model errors below the topography and in its top 50 m, and traveltime residuals
top = np.argmax(mask, axis=0)
shallow = mask & (np.arange(nz)[:, None] - top[None, :] < 5)
t = traveltimes('data')
print(f'\n{"Model":<30}{"RMS error (m/s)":>17}{"Top 50 m (m/s)":>16}{"Residual (ms)":>15}')
rows = [('Initial model', vp0, 'test/iteration_0/synthetic')] + \
    [(label, model(run, niter), f'{run}/iteration_{niter}/synthetic') for run, label, *_ in runs]
for label, m, synthetic in rows:
    print(f'{label:<30}{rms((m - vp)[mask]):>17.1f}{rms((m - vp)[shallow]):>16.1f}'
          f'{rms(traveltimes(synthetic) - t)*1000:>15.2f}')

Path('figures').mkdir(exist_ok=True)

# Models; the region above the topography is blank
extent = [-0.5*d, (nx - 0.5)*d, (nz - 0.5)*d, -0.5*d]
panels = [('True model', vp), ('Initial model', vp0)] + [(label, model(run, niter)) for run, label, *_ in runs[:4]]
fig, axes = plt.subplots(3, 2, figsize=(10.0, 5.6), constrained_layout=True)
for ax, (title, m), tag in zip(axes.ravel(), panels, 'abcdef'):
    im = ax.imshow(np.where(mask, m, np.nan), extent=extent, cmap='jet', vmin=2000, vmax=6000,
                   interpolation='none', aspect='auto')
    ax.set_title(f'({tag}) {title}', loc='left')
    ax.set_yticks([0, 0.5, 1])
for ax in axes[:, 0]:
    ax.set_ylabel('Depth (km)')
for ax in axes[:, 1]:
    ax.tick_params(labelleft=False)
for ax in axes[:2].ravel():
    ax.tick_params(labelbottom=False)
for ax in axes[2]:
    ax.set_xlabel('Distance (km)')
cb = fig.colorbar(im, ax=axes, location='right', shrink=0.8, aspect=30)
cb.set_label('P-wave velocity (m/s)')
fig.savefig('figures/regularization_models.pdf')
fig.savefig('figures/regularization_models.png', dpi=200)

# Data misfit, normalized by its value in the initial model, and model error in each iteration
fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.6), constrained_layout=True)
for run, label, _, color, style in runs:
    misfit = np.loadtxt(f'{run}/data_misfit.txt')
    axes[0].semilogy(misfit[:, 0], misfit[:, 2], color=color, ls=style, lw=1.8, label=label)
    axes[1].plot(range(niter + 1), [rms((model(run, k) - vp)[mask]) for k in range(niter + 1)],
                 color=color, ls=style, lw=1.8, label=label)
axes[0].set_ylabel('Normalized data misfit')
axes[0].set_title('(a) Data misfit', loc='left')
axes[0].grid(True, which='both', color='0.88', lw=0.6)
axes[1].set_ylabel('RMS model error (m/s)')
axes[1].set_title('(b) Model error', loc='left')
axes[1].grid(True, color='0.88', lw=0.6)
for ax in axes:
    ax.set_xlabel('Iteration')
    ax.set_xlim(0, niter)
axes[1].legend(frameon=False, loc='upper right')
fig.savefig('figures/regularization_curves.pdf')
fig.savefig('figures/regularization_curves.png', dpi=200)
