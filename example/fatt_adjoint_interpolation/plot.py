"""Compare the gradients of both adjoint_receiver settings with the
finite-difference reference, print the errors and plot the figures."""

import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

plt.rcParams.update({"font.family": "Arial", "font.size": 10, "savefig.dpi": 300, "savefig.bbox": "tight",
                     "savefig.pad_inches": 0.02})

setup = json.load(open("setup.json"))
nx, nz, d = setup["nx"], setup["nz"], setup["d"]
receivers = setup["receivers"]
windows = setup["windows"]
v = np.fromfile("model/v.bin", dtype=np.float32).reshape(nx, nz).T
nodes = np.load("results/fd_nodes.npy")
fd = np.load("results/fd_kernels.npy")
extent = [-d/2, (nx - 0.5)*d, (nz - 0.5)*d, -d/2]


def in_window(name):
    w = windows[name]
    return (nodes[:, 0] >= w[0]) & (nodes[:, 0] <= w[1]) & (nodes[:, 1] >= w[2]) & (nodes[:, 1] <= w[3])


def adjoint(case, mode):
    """The adjoint gradient of 1/2 sum e^2 with respect to the velocity at each grid point,
    -lambda/v^3 times the cell area, at the window points."""
    g = np.load(f"results/grad_{case}_{mode}.npy")*d*d
    return g[nodes[:, 1], nodes[:, 0]]


def compare(g, k):
    """The ratio of the adjoint to the reference amplitude, and the relative errors
    of the adjoint gradient as is and after scaling it to the reference."""
    ratio = (g @ k)/(k @ k)
    scale = (g @ k)/(g @ g) if g @ g > 0 else 0.0
    return ratio, np.linalg.norm(g - k)/np.linalg.norm(k), np.linalg.norm(scale*g - k)/np.linalg.norm(k)


def image(ax, values, mask, window, **kwargs):
    a = np.full((nz, nx), np.nan)
    a[nodes[mask, 1], nodes[mask, 0]] = values
    im = ax.imshow(a, extent=extent, aspect="equal", interpolation="nearest", **kwargs)
    ax.set_xlim(window[0]*d - d/2, window[1]*d + d/2)
    ax.set_ylim(window[3]*d + d/2, window[2]*d - d/2)
    return im


def label(ax, s):
    ax.text(0.03, 0.96, s, transform=ax.transAxes, ha="left", va="top", fontweight="bold",
            bbox=dict(facecolor="white", edgecolor="none", pad=0.8, alpha=0.8))


# ---------------------------------------------------------------- errors
sweep = [k for k, r in enumerate(receivers) if r["case"].startswith("sweep")]
pair = [k for k, r in enumerate(receivers) if r["case"] == "pair"]
ms, mp = in_window("sweep"), in_window("pair")
rows = []
print("moving receiver: shape error of the kernel near the receiver against finite differences")
print(" offset   linear  nearest")
for k in sweep:
    case = receivers[k]["case"]
    ref = fd[k, ms]
    lin = compare(adjoint(case, "linear")[ms]/receivers[k]["residual"], ref)
    near = compare(adjoint(case, "nearest")[ms]/receivers[k]["residual"], ref)
    rows.append((receivers[k]["z"] - 200.0, lin, near))
    print(f"  {receivers[k]['z'] - 200.0:4.0f} m  {100*lin[2]:6.1f}%  {100*near[2]:6.1f}%")
ref_pair = np.array([receivers[k]["residual"] for k in pair]) @ fd[pair][:, mp]
print("receiver pair (+1 ms, -1 ms in one cell): gradient against finite differences")
for mode in ["linear", "nearest"]:
    g = adjoint("pair", mode)[mp]
    if np.abs(g).max() == 0:
        print(f"  {mode:8s} zero gradient")
    else:
        r = compare(g, ref_pair)
        print(f"  {mode:8s} shape error {100*r[2]:.1f}%")
json.dump({"sweep": [[a, list(b), list(c)] for a, b, c in rows]}, open("results/errors.json", "w"), indent=1)

# ---------------------------------------------------------------- setup
fig, ax = plt.subplots(figsize=(6.5, 3.2))
im = ax.imshow(v, extent=extent, cmap="jet", aspect="equal")
ax.plot(*setup["source"], "*", color="red", mec="black", ms=14, clip_on=False, zorder=5)
for name, w in windows.items():
    ax.add_patch(Rectangle((w[0]*d - d/2, w[2]*d - d/2), (w[1] - w[0] + 1)*d, (w[3] - w[2] + 1)*d,
                           fill=False, edgecolor="black", lw=1.2))
ax.plot([receivers[k]["x"] for k in sweep], [receivers[k]["z"] for k in sweep], "o", color="white", mec="black", ms=3)
ax.plot([receivers[k]["x"] for k in pair], [receivers[k]["z"] for k in pair], "o", color="white", mec="black", ms=3)
ax.text(1015, 170, "moving receiver", ha="right", va="bottom")
ax.text(835, 365, "receiver pair", ha="right", va="bottom")
ax.set_xlabel("Horizontal position (m)")
ax.set_ylabel("Depth (m)")
cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.015)
cb.set_label("Velocity (m/s)")
fig.savefig("fig_setup.png")
plt.close(fig)

# ---------------------------------------------------------------- kernels and pair gradients
k5 = sweep[5]
case5 = receivers[k5]["case"]
fig, axes = plt.subplots(2, 3, figsize=(6.5, 4.4), gridspec_kw=dict(wspace=0.08, hspace=0.32))
ref = -fd[k5, ms]*1e6
lin = -adjoint(case5, "linear")[ms]/1e-3*1e6
near = -adjoint(case5, "nearest")[ms]/1e-3*1e6
near *= (near @ ref)/(near @ near)
vmax = ref.max()
for c, (a, title) in enumerate([(ref, "Finite difference"), (lin, "Linear"), (near, "Nearest (scaled)")]):
    im = image(axes[0, c], a, ms, windows["sweep"], cmap="hot_r", vmin=0, vmax=vmax)
    axes[0, c].set_title(title, fontsize=10)
    axes[0, c].plot(receivers[k5]["x"], receivers[k5]["z"], "x", color="blue", ms=9, mew=2)
pos = axes[0, 2].get_position()
cb = fig.colorbar(im, cax=fig.add_axes([pos.x1 + 0.012, pos.y0, 0.012, pos.height]))
cb.set_label(r"$-\partial T/\partial v$ ($\mu$s per m/s)", fontsize=9)
axes[0, 2].plot(1000.0, 210.0, "o", mfc="none", mec="blue", ms=10, mew=1.5)
gl = adjoint("pair", "linear")[mp]
gn = adjoint("pair", "nearest")[mp]
vmax = np.abs(ref_pair).max()*1e10
for c, a in enumerate([ref_pair, gl, gn]):
    im = image(axes[1, c], a*1e10, mp, windows["pair"], cmap="seismic", vmin=-vmax, vmax=vmax)
    axes[1, c].plot(receivers[pair[0]]["x"], receivers[pair[0]]["z"], "+", color="black", ms=9, mew=2)
    axes[1, c].plot(receivers[pair[1]]["x"], receivers[pair[1]]["z"], "_", color="black", ms=9, mew=2)
if np.abs(gn).max() == 0:
    axes[1, 2].text(0.5, 0.25, "zero", transform=axes[1, 2].transAxes, ha="center", va="center", fontsize=12)
pos = axes[1, 2].get_position()
cb = fig.colorbar(im, cax=fig.add_axes([pos.x1 + 0.012, pos.y0, 0.012, pos.height]))
cb.set_label(r"$\partial J/\partial v$ ($10^{-10}$ s$^2$ per m/s)", fontsize=9)
for r in range(2):
    axes[r, 0].set_ylabel("Depth (m)")
    for c in range(3):
        axes[r, c].set_xlabel("x (m)")
        if c > 0:
            axes[r, c].set_yticklabels([])
        label(axes[r, c], "(" + "abcdef"[3*r + c] + ")")
fig.savefig("fig_gradients.png")
plt.close(fig)

# ---------------------------------------------------------------- moving receiver
offset = np.array([r[0] for r in rows])
fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.6), gridspec_kw=dict(wspace=0.35))
axes[0].plot(offset, [100*r[2][2] for r in rows], "o-", color="red", ms=4, label="Nearest")
axes[0].plot(offset, [100*r[1][2] for r in rows], "s-", color="blue", ms=4, label="Linear")
axes[0].set_ylabel("Shape error (%)")
axes[0].set_ylim(0, 1.35*max(100*r[2][2] for r in rows))
axes[0].legend(loc="upper right", fontsize=9, frameon=False)
near_amp = np.array([r[2][0] for r in rows])
axes[1].plot(offset, near_amp/near_amp[0], "o-", color="red", ms=4, label="Nearest, relative to 0 m")
axes[1].plot(offset, [r[1][0] for r in rows], "s-", color="blue", ms=4, label="Linear")
axes[1].axhline(1, color="black", lw=0.6, ls=":")
amp = np.concatenate([near_amp/near_amp[0], [r[1][0] for r in rows]])
axes[1].set_ylim(amp.min() - 0.55*np.ptp(amp), amp.max() + 0.08*np.ptp(amp))
axes[1].set_ylabel("Relative amplitude")
axes[1].legend(loc="lower right", fontsize=9, frameon=False)
for c in range(2):
    axes[c].set_xlabel("Receiver offset (m)")
    axes[c].set_xticks([0, 2, 4, 6, 8, 10])
    axes[c].axvline(5, color="black", lw=0.6, ls=":")
    label(axes[c], "(" + "ab"[c] + ")")
fig.savefig("fig_moving_receiver.png")
plt.close(fig)
print("figures: fig_setup.png, fig_gradients.png, fig_moving_receiver.png")
