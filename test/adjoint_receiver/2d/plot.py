"""Figures of the 2D kernel test."""

import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyze import meta, read_array, grad, load_stack, nx, nz, d, v, X, Z

plt.rcParams.update({"font.family": "Arial", "font.size": 10, "axes.linewidth": 0.8,
                     "xtick.major.width": 0.8, "ytick.major.width": 0.8,
                     "savefig.dpi": 300, "savefig.bbox": "tight", "savefig.pad_inches": 0.02})
idx = meta["index"]
metr = json.load(open("out/metrics.json"))
extent = [-d/2, (nx - 0.5)*d, (nz - 0.5)*d, -d/2]


def label(ax, s):
    ax.text(0.02, 0.97, s, transform=ax.transAxes, ha="left", va="top", fontsize=10, fontweight="bold",
            bbox=dict(facecolor="white", edgecolor="none", pad=1.0, alpha=0.85))


# ---------------------------------------------------------------- setup
fig, ax = plt.subplots(figsize=(6.5, 3.1))
im = ax.imshow(v, extent=extent, cmap="jet", aspect="equal")
g = meta["groups"]
ax.plot(*meta["src1"], marker="*", color="red", ms=14, mec="black", clip_on=False, zorder=5)
ax.plot(*meta["src2"], marker="*", color="yellow", ms=14, mec="black", clip_on=False, zorder=5)
ax.plot(*np.array(g["surface"]).T, "v", color="white", mec="black", ms=5, clip_on=False, zorder=4)
for nm in ["sweep_x", "sweep_z", "pair_same", "pair_adj"]:
    p = np.array(g[nm])
    ax.plot(p[:, 0], p[:, 1], "o", color="white", mec="black", ms=3.5)
ax.plot(*np.array(meta["vsp"]).T, "s", color="yellow", mec="black", ms=2.5)
ax.annotate("sweep along x", (910, 300), (1010, 380), ha="center", arrowprops=dict(arrowstyle="-", color="black"))
ax.annotate("sweep along z", (1000, 205), (1060, 130), ha="center", arrowprops=dict(arrowstyle="-", color="black"))
ax.annotate("receiver pairs", (803, 403), (830, 500), ha="left", arrowprops=dict(arrowstyle="-", color="black"))
ax.text(712, 560, "borehole", ha="left", va="bottom")
ax.set_xlabel("Horizontal position (m)")
ax.set_ylabel("Depth (m)")
cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.015)
cb.set_label("Velocity (m/s)")
fig.savefig("fig_setup.png")
plt.close(fig)

# ---------------------------------------------------------------- kernels near a receiver
trec1, fd1, near1, lin1 = (np.fromfile("out/fd_shot1.bin.trec", dtype=np.float32),
                           *(load_stack(f, 90) for f in ["out/fd_shot1.bin", "out/each_shot1_nearest.bin",
                                                          "out/each_shot1_linear.bin"]))
cases = [("sweep_x", 5), ("sweep_z", 5)]
fig, axes = plt.subplots(2, 3, figsize=(6.5, 4.55), sharex="row", sharey="row",
                         gridspec_kw=dict(wspace=0.06, hspace=0.32))
for row, (nm, s) in enumerate(cases):
    r = idx[nm][s]
    rx, rz = meta["groups"][nm][s]
    k = -fd1[r]*1e6
    gl = -grad(lin1[r])*1e6
    alpha = metr[f"shot1_{r}"]["near_all"][1]
    gn = -grad(near1[r])*alpha*1e6
    vmax = k[np.hypot(X - rx, Z - rz) < 60].max()
    ix0, iz0 = int(np.floor(rx/d + 0.5)), int(np.floor(rz/d + 0.5))
    for col, (a, ttl) in enumerate([(k, "Finite difference"), (gl, "Linear"), (gn, "Nearest (rescaled)")]):
        ax = axes[row, col]
        im = ax.imshow(a, extent=extent, cmap="hot_r", vmin=0, vmax=vmax, aspect="equal", interpolation="nearest")
        ax.set_xlim(rx - 62, rx + 18)
        ax.set_ylim(rz + 32, rz - 48)
        ax.set_xticks(np.arange(np.ceil((rx - 62)/20), np.floor((rx + 18)/20) + 1)*20)
        gx = np.arange(np.ceil((rx - 62)/d), np.floor((rx + 18)/d) + 1)*d
        gz = np.arange(np.ceil((rz - 48)/d), np.floor((rz + 32)/d) + 1)*d
        GX, GZ = np.meshgrid(gx, gz)
        ax.plot(GX.ravel(), GZ.ravel(), ".", color="black", ms=1.5)
        ax.plot(rx, rz, "x", color="blue", ms=9, mew=2)
        if col == 2:
            ax.plot(ix0*d, iz0*d, "o", mfc="none", mec="blue", ms=10, mew=1.5)
        if row == 0:
            ax.set_title(ttl, fontsize=10)
        if col == 0:
            ax.set_ylabel("Depth (m)")
        ax.set_xlabel("Horizontal position (m)")
        label(ax, "(" + "abcdef"[3*row + col] + ")")
    cax = fig.add_axes([axes[row, 2].get_position().x1 + 0.012, axes[row, 2].get_position().y0,
                        0.012, axes[row, 2].get_position().height])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(r"$-\partial T/\partial v$ ($\mu$s per m/s)")
fig.savefig("fig_kernels.png")
plt.close(fig)

# ---------------------------------------------------------------- sub-grid sweeps
fig, axes = plt.subplots(2, 2, figsize=(6.5, 4.4), sharex=True, gridspec_kw=dict(wspace=0.32, hspace=0.12))
for row, nm in enumerate(["sweep_x", "sweep_z"]):
    rows = [metr[f"shot1_{r}"] for r in idx[nm]]
    s = np.arange(len(rows))
    en = np.array([r["near_rec"][3] for r in rows])*100
    el = np.array([r["lin_rec"][2] for r in rows])*100
    els = np.array([r["lin_rec"][3] for r in rows])*100
    an = np.array([r["near_all"][1] for r in rows])
    al = np.array([r["lin_all"][1] for r in rows])
    ax = axes[row, 0]
    ax.plot(s, en, "o-", color="red", ms=4, label="Nearest (rescaled)")
    ax.plot(s, el, "s-", color="blue", ms=4, label="Linear")
    ax.plot(s, els, "^--", color="blue", ms=4, mfc="none", label="Linear (rescaled)")
    ax.set_ylabel("Error near receiver (%)")
    ax.set_ylim(0, None)
    ax = axes[row, 1]
    ax.plot(s, an/an[0], "o-", color="red", ms=4, label="Nearest, relative to its value at 0 m")
    ax.plot(s, al, "s-", color="blue", ms=4, label="Linear")
    ax.set_ylabel("Least-squares factor")
    ax.axhline(1, color="black", lw=0.6, ls=":")
    for c in range(2):
        axes[row, c].axvline(5, color="black", lw=0.6, ls=":")
        axes[row, c].axvline(15, color="black", lw=0.6, ls=":")
        label(axes[row, c], "(" + "abcd"[2*row + c] + ")")
    axes[row, 0].text(0.98, 0.95, "along x" if row == 0 else "along z", transform=axes[row, 0].transAxes,
                      ha="right", va="top")
axes[0, 0].legend(loc="lower left", fontsize=9, frameon=False, bbox_to_anchor=(0.0, 1.0), ncol=2,
                  columnspacing=0.8, handlelength=1.6)
axes[0, 1].legend(loc="lower left", fontsize=9, frameon=False, bbox_to_anchor=(0.0, 1.0), ncol=1,
                  handlelength=1.6)
for c in range(2):
    axes[1, c].set_xlabel("Receiver offset from a grid point (m)")
    axes[1, c].set_xticks([0, 5, 10, 15, 20])
fig.savefig("fig_sweep.png")
plt.close(fig)

# ---------------------------------------------------------------- group gradients
gm = json.load(open("out/group_metrics.json"))


def group_row(axes, gname, window, letters, scale=1e10):
    gfd = np.load(f"out/gfd_{gname}.npy")
    gl = grad(read_array(f"out/adj_{gname}_linear.bin", (nz, nx)).astype(float))
    gn = grad(read_array(f"out/adj_{gname}_nearest.bin", (nz, nx)).astype(float))
    alpha = gm[f"{gname}_nearest"]["alpha"]
    gn = gn*(alpha if np.isfinite(alpha) else 1.0)
    x0, x1, z0, z1 = window
    inwin = (X >= x0) & (X <= x1) & (Z >= z0) & (Z <= z1)
    src = meta["src2"] if gname == "vsp" else meta["src1"]
    inwin &= np.hypot(X - src[0], Z - src[1]) > 50
    vmax = np.percentile(np.abs(gfd[inwin]), 99.5)*scale
    for col, a in enumerate([gfd, gl, gn]):
        ax = axes[col]
        im = ax.imshow(a*scale, extent=extent, cmap="seismic", vmin=-vmax, vmax=vmax, aspect="equal",
                       interpolation="nearest")
        ax.set_xlim(x0, x1)
        ax.set_ylim(z1, z0)
        label(ax, letters[col])
    return im


fig, axes = plt.subplots(2, 3, figsize=(6.5, 4.6), sharex="row", sharey="row",
                         gridspec_kw=dict(wspace=0.06, hspace=0.3))
for row, gname in enumerate(["pair_same", "pair_adj"]):
    p = np.array(meta["groups"][gname])
    c = p.mean(axis=0)
    im = group_row(axes[row], gname, (c[0] - 60, c[0] + 40, c[1] - 70, c[1] + 30),
                   ["(" + q + ")" for q in "abcdef"[3*row:3*row + 3]], scale=1e10)
    for col in range(3):
        ax = axes[row, col]
        gx = np.arange(np.ceil((c[0] - 60)/d), np.floor((c[0] + 40)/d) + 1)*d
        gz = np.arange(np.ceil((c[1] - 70)/d), np.floor((c[1] + 30)/d) + 1)*d
        GX, GZ = np.meshgrid(gx, gz)
        ax.plot(GX.ravel(), GZ.ravel(), ".", color="black", ms=1.5)
        ax.plot(p[0, 0], p[0, 1], "+", color="black", ms=9, mew=2)
        ax.plot(p[1, 0], p[1, 1], "_", color="black", ms=9, mew=2)
        ax.set_xlabel("Horizontal position (m)")
        if row == 0:
            ax.set_title(["Finite difference", "Linear", "Nearest (rescaled)"][col], fontsize=10)
    axes[row, 0].set_ylabel("Depth (m)")
    pos = axes[row, 2].get_position()
    cb = fig.colorbar(im, cax=fig.add_axes([pos.x1 + 0.012, pos.y0, 0.012, pos.height]))
    cb.set_label(r"$\partial J/\partial v$ ($10^{-10}$ s$^2$ per m/s)", fontsize=9)
axes[0, 2].text(0.5, 0.5, "zero", transform=axes[0, 2].transAxes, ha="center", va="center", fontsize=12)
fig.savefig("fig_pairs.png")
plt.close(fig)

def group_fields(gname):
    gfd = np.load(f"out/gfd_{gname}.npy")
    gl = grad(read_array(f"out/adj_{gname}_linear.bin", (nz, nx)).astype(float))
    gn = grad(read_array(f"out/adj_{gname}_nearest.bin", (nz, nx)).astype(float))*gm[f"{gname}_nearest"]["alpha"]
    return gfd, gl, gn


# borehole: true gradient, errors of both modes, and profiles
gfd, gl, gn = group_fields("vsp")
scale = 1e7
x0, x1, z0, z1 = 600, 730, 0, 560
win = (X >= x0) & (X <= x1) & (np.hypot(X - meta["src2"][0], Z - meta["src2"][1]) > 50)
vmax = np.percentile(np.abs(gfd[win]), 99.5)*scale
fig = plt.figure(figsize=(6.5, 4.0))
gs = fig.add_gridspec(1, 6, width_ratios=[1, 1, 1, 0.07, 0.55, 1.45], wspace=0.08)
axes = [fig.add_subplot(gs[0, c]) for c in range(3)]
for c, (a, ttl) in enumerate([(gfd, "Finite difference"), (gl - gfd, "Linear $-$ FD"),
                              (gn - gfd, "Nearest $-$ FD")]):
    ax = axes[c]
    im = ax.imshow(a*scale, extent=extent, cmap="seismic", vmin=-vmax, vmax=vmax, aspect="auto",
                   interpolation="nearest")
    ax.plot(*np.array(meta["vsp"]).T, ".", color="black", ms=1.2)
    ax.set_xlim(x0, x1)
    ax.set_ylim(z1, z0)
    ax.set_xticks([620, 700])
    ax.set_title(ttl, fontsize=10)
    ax.set_xlabel("x (m)")
    if c > 0:
        ax.set_yticklabels([])
    label(ax, "(" + "abc"[c] + ")")
axes[0].set_ylabel("Depth (m)")
cb = fig.colorbar(im, cax=fig.add_subplot(gs[0, 3]))
cb.set_label(r"$\partial J/\partial v$ ($10^{-7}$ s$^2$ per m/s)", fontsize=9)
ax = fig.add_subplot(gs[0, 5])
iz = int(300/d)
xg = np.arange(nx)*d
ax.plot(xg, gfd[iz]*scale, "-", color="black", lw=2.2, label="Finite difference")
ax.plot(xg, gl[iz]*scale, "--", color="blue", lw=1.5, label="Linear")
ax.plot(xg, gn[iz]*scale, "-", color="red", lw=1.2, label="Nearest")
ax.set_xlim(640, 720)
ax.set_xlabel("x (m)")
ax.set_ylabel(r"$\partial J/\partial v$ at 300 m depth ($10^{-7}$)", fontsize=9)
ax.legend(loc="upper left", fontsize=9, frameon=False, bbox_to_anchor=(0.0, 0.92), handlelength=1.5)
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")
label(ax, "(d)")
fig.savefig("fig_vsp.png")
plt.close(fig)

# surface line: true gradient and errors of both modes
gfd, gl, gn = group_fields("surface")
scale = 1e8
win = (Z <= 150) & (np.hypot(X - meta["src1"][0], Z - meta["src1"][1]) > 50)
vmax = np.percentile(np.abs(gfd[win]), 99.5)*scale
fig, axes = plt.subplots(3, 1, figsize=(6.5, 3.9), sharex=True, gridspec_kw=dict(hspace=0.12))
for r, (a, ttl) in enumerate([(gfd, "Finite difference"), (gl - gfd, "Linear $-$ FD"),
                              (gn - gfd, "Nearest $-$ FD")]):
    ax = axes[r]
    im = ax.imshow(a*scale, extent=extent, cmap="seismic", vmin=-vmax, vmax=vmax, aspect="auto",
                   interpolation="nearest")
    ax.plot(*np.array(meta["groups"]["surface"]).T, "v", color="black", ms=2.5, clip_on=False)
    ax.set_ylim(150, 0)
    ax.set_ylabel("Depth (m)")
    ax.text(0.99, 0.08, ttl, transform=ax.transAxes, ha="right", va="bottom",
            bbox=dict(facecolor="white", edgecolor="none", pad=1.0, alpha=0.85))
    label(ax, "(" + "abc"[r] + ")")
axes[-1].set_xlabel("Horizontal position (m)")
pos0, pos2 = axes[0].get_position(), axes[2].get_position()
cb = fig.colorbar(im, cax=fig.add_axes([pos0.x1 + 0.012, pos2.y0, 0.015, pos0.y1 - pos2.y0]))
cb.set_label(r"$\partial J/\partial v$ ($10^{-8}$ s$^2$ per m/s)")
fig.savefig("fig_surface.png")
plt.close(fig)

print("done")
