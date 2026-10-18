"""Errors and figures of the 2D TRTT example."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({"font.family": "Arial", "font.size": 10, "savefig.dpi": 300, "savefig.bbox": "tight",
                     "savefig.pad_inches": 0.02})
nx, nz, d = 301, 151, 10.0
x, z = np.arange(nx)*d/1e3, np.arange(nz)*d/1e3
ext = [0, x[-1], z[-1], 0]
lit = slice(15, 286)                    # columns 150 m from the ends
runs = [("first_arrivals", "First arrivals only", "blue", "-."),
        ("trtt", "TRTT", "red", "-"),
        ("true_reflectors", "True reflector fixed", "black", ":")]


def load(path, shape=(nz, nx)):
    return np.fromfile(path, dtype=np.float32).reshape(shape[::-1]).T


def iterations(run):
    return sorted(int(p.name.split("_")[1]) for p in Path(f"test_{run}").glob("iteration_*")
                  if (p/"model"/"updated_vp.bin").is_file())


def depth_of(refl, label=1):
    """The depth of a reflector in each column of a label-plus-offset model."""
    iz = np.floor(refl + 0.5)
    zz = np.where(iz == label, (np.arange(nz)[:, None] + refl - label)*d, np.nan)
    with np.errstate(all="ignore"):
        return np.nanmean(zz, axis=0)


vp_true, vp_init, depth = load("model/vp_true.bin"), load("model/vp_init.bin"), np.load("model/depth.npy")
rms = lambda w: float(np.sqrt(np.mean(w**2)))
results = {}
for run, name, color, ls in runs:
    its = iterations(run)
    vp = [vp_init] + [load(f"test_{run}/iteration_{k}/model/updated_vp.bin") for k in its]
    out = {"vp": vp[-1], "error": [rms(v - vp_true) for v in vp],
           "misfit": np.loadtxt(f"test_{run}/data_misfit.txt", ndmin=2)[:, 2], "last": its[-1]}
    if run == "trtt":
        refl = [load(f"test_{run}/iteration_{k}/model/updated_refl.bin") for k in its]
        out["depth"] = depth_of(refl[-1])
        out["refl_error"] = [float(np.nanmean(np.abs(depth_of(r)[lit] - depth[lit]))) for r in refl]
    results[run] = out

# Velocity change from the starting model, with the reflector
fig, axes = plt.subplots(2, 2, figsize=(10, 5.0), sharex=True, sharey=True, gridspec_kw=dict(hspace=0.25, wspace=0.05))
panels = [("True model", vp_true, None)] + [(name, results[run]["vp"], run) for run, name, *_ in runs]
panels = [panels[0], panels[2], panels[1], panels[3]]
for ax, (name, v, run), s in zip(axes.flat, panels, "abcd"):
    im = ax.imshow(v - vp_init, extent=ext, cmap="seismic", vmin=-300, vmax=300, aspect="auto", interpolation="none")
    ax.plot(x, depth/1e3, "--" if run == "trtt" else "-", color="black", lw=1.0)
    if run == "trtt":
        ax.plot(x, results[run]["depth"]/1e3, "-", color="black", lw=1.2)
    title = name if run is None else f"{name}: {results[run]['error'][-1]:.1f} m/s"
    ax.set_title(title, fontsize=10)
    ax.text(0.015, 0.95, f"({s})", transform=ax.transAxes, ha="left", va="top", fontweight="bold")
for ax in axes[:, 0]:
    ax.set_ylabel("z (km)")
for ax in axes[-1, :]:
    ax.set_xlabel("x (km)")
    ax.set_xticks([0, 0.5, 1, 1.5, 2, 2.5])
cb = fig.colorbar(im, ax=axes, fraction=0.03, pad=0.015)
cb.set_label("Vp change (m/s)")
fig.savefig("fig_models.png")
plt.close(fig)

# Convergence
fig, axes = plt.subplots(1, 3, figsize=(10, 3.2), gridspec_kw=dict(wspace=0.32))
for run, name, color, ls in runs:
    o = results[run]
    axes[0].plot(np.arange(len(o["error"])), o["error"], ls, color=color, lw=1.5, label=name)
    axes[2].semilogy(np.arange(len(o["misfit"])), o["misfit"], ls, color=color, lw=1.5)
axes[1].plot(np.arange(1, len(results["trtt"]["refl_error"]) + 1), results["trtt"]["refl_error"], "-", color="red", lw=1.5)
for ax, label, s in zip(axes, ["Vp RMS error (m/s)", "Reflector |error| (m)", "Normalized misfit"], "abc"):
    ax.set_xlabel("Iteration")
    ax.set_ylabel(label)
    ax.text(0.03, 0.95, f"({s})", transform=ax.transAxes, ha="left", va="top", fontweight="bold")
axes[1].set_ylim(0, None)
fig.legend(*axes[0].get_legend_handles_labels(), loc="lower center", bbox_to_anchor=(0.5, 0.95), ncol=3, frameon=False)
fig.savefig("fig_convergence.png")
plt.close(fig)

# Traveltime residuals of three shots, after the last iteration
xr = np.arange(0.0, 3001.0, 20.0)/1e3
shots = {3: 0.5, 8: 1.5, 13: 2.5}
fig, axes = plt.subplots(2, 3, figsize=(10, 5.0), sharex=True, gridspec_kw=dict(hspace=0.12, wspace=0.25))
for col, (sid, sx) in enumerate(shots.items()):
    obs = np.fromfile(f"data/shot_{sid}_traveltime_p.bin", dtype=np.float32).reshape(2, -1).T
    for run, name, color, ls in runs:
        syn = np.fromfile(f"test_{run}/iteration_{results[run]['last']}/synthetic/shot_{sid}_traveltime_p.bin",
                          dtype=np.float32).reshape(2, -1).T
        for e in range(2):
            axes[e, col].plot(xr, (syn[:, e] - obs[:, e])*1e3, ls, color=color, lw=1.2, label=name)
    axes[0, col].set_title(f"Source at x = {sx:g} km", fontsize=10)
    for e in range(2):
        axes[e, col].axhline(0, color="black", lw=0.6)
        axes[e, col].text(0.03, 0.95, f"({'abcdef'[3*e + col]})", transform=axes[e, col].transAxes, ha="left",
                          va="top", fontweight="bold")
    axes[1, col].set_xlabel("Receiver x (km)")
axes[0, 0].set_ylabel("First arrival\nresidual (ms)")
axes[1, 0].set_ylabel("Reflection\nresidual (ms)")
fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="lower center", bbox_to_anchor=(0.5, 0.93), ncol=3,
           frameon=False)
fig.savefig("fig_traveltimes.png")
plt.close(fig)

print(f"Vp RMS error, from {results['trtt']['error'][0]:.1f} m/s at the start:")
for run, name, *_ in runs:
    print(f"  {name:22s} {results[run]['error'][-1]:5.1f} m/s, normalized misfit {results[run]['misfit'][-1]:.2e}")
print(f"TRTT reflector mean |error|: {results['trtt']['refl_error'][-1]:.1f} m")
print("Figures: fig_models.png, fig_convergence.png, fig_traveltimes.png")
