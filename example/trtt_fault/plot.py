"""Errors and figures of the 2D TRTT example with faults."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({"font.family": "Arial", "font.size": 10, "savefig.dpi": 300, "savefig.bbox": "tight",
                     "savefig.pad_inches": 0.02})
nx, nz, d, nrefl = 401, 161, 10.0, 3
x, z = np.arange(nx)*d/1e3, np.arange(nz)*d/1e3
ext = [0, x[-1], z[-1], 0]
lit = slice(20, 381)                    # columns x = 0.2 to 3.8 km
runs = [("first_arrivals", "First arrivals only", "blue", "-."),
        ("trtt_all_offsets", "TRTT, reflectors from all offsets", "magenta", "--"),
        ("trtt", "TRTT, reflectors from offsets up to 1 km", "red", "-"),
        ("true_reflectors", "True reflectors fixed", "black", ":")]


def load(path, shape=(nz, nx)):
    return np.fromfile(path, dtype=np.float32).reshape(shape[::-1]).T


def iterations(run):
    return sorted(int(p.name.split("_")[1]) for p in Path(f"test_{run}").glob("iteration_*")
                  if (p/"model"/"updated_vp.bin").is_file())


def depth_of(refl, label):
    """The depth of a reflector in each column of a label-plus-offset model."""
    iz = np.floor(refl + 0.5)
    zz = np.where(iz == label, (np.arange(nz)[:, None] + refl - label)*d, np.nan)
    with np.errstate(all="ignore"):
        return np.nanmean(zz, axis=0)


vp_true, vp_init, depth = load("model/vp_true.bin"), load("model/vp_init.bin"), np.load("model/depth.npy")
# The target: the grid points above the deepest reflector in the lit columns
target = np.zeros((nz, nx), dtype=bool)
for j in range(lit.start, lit.stop):
    target[np.arange(nz)*d <= depth[-1][j], j] = True
rms = lambda w: float(np.sqrt(np.mean(w**2)))
results = {}
for run, name, color, ls in runs:
    its = iterations(run)
    vp = [vp_init] + [load(f"test_{run}/iteration_{k}/model/updated_vp.bin") for k in its]
    out = {"vp": vp[-1], "error": [rms((v - vp_true)[target]) for v in vp],
           "misfit": np.loadtxt(f"test_{run}/data_misfit.txt", ndmin=2)[:, 2], "last": its[-1]}
    if run.startswith("trtt"):
        refl = [load(f"test_{run}/iteration_{k}/model/updated_refl.bin") for k in its]
        out["depth"] = [depth_of(refl[-1], l) for l in range(1, nrefl + 1)]
        out["refl_error"] = np.array([[float(np.nanmean(np.abs(depth_of(r, l)[lit] - depth[l - 1][lit])))
                                       for l in range(1, nrefl + 1)] for r in refl])
    results[run] = out

# The velocity models with the reflectors
fig, axes = plt.subplots(3, 2, figsize=(10, 6.4), sharex=True, sharey=True, gridspec_kw=dict(hspace=0.3, wspace=0.08))
title = {run: f"{name}: {results[run]['error'][-1]:.1f} m/s" for run, name, *_ in runs}
# Each panel: title, velocity, and the reflectors drawn on it: the true ones,
# or those of a TRTT run (solid) over the true ones (dashed)
panels = [("True model", vp_true, "true"),
          (f"Starting model: {results['trtt']['error'][0]:.1f} m/s", vp_init, None),
          (title["first_arrivals"], results["first_arrivals"]["vp"], None),
          (title["true_reflectors"], results["true_reflectors"]["vp"], "true"),
          (title["trtt_all_offsets"], results["trtt_all_offsets"]["vp"], "trtt_all_offsets"),
          (title["trtt"], results["trtt"]["vp"], "trtt")]
for ax, (name, v, refl), s in zip(axes.flat, panels, "abcdef"):
    im = ax.imshow(v, extent=ext, cmap="jet", vmin=1800, vmax=3600, aspect="auto", interpolation="none")
    for l in range(nrefl):
        if refl == "true":
            ax.plot(x, depth[l]/1e3, "-", color="black", lw=1.0)
        elif refl is not None:
            ax.plot(x, depth[l]/1e3, "--", color="black", lw=0.8)
            ax.plot(x, results[refl]["depth"][l]/1e3, "-", color="black", lw=1.2)
    ax.set_title(name, fontsize=10)
    ax.text(0.012, 0.93, f"({s})", transform=ax.transAxes, ha="left", va="top", fontweight="bold", color="white")
for ax in axes[:, 0]:
    ax.set_ylabel("z (km)")
for ax in axes[-1, :]:
    ax.set_xlabel("x (km)")
    ax.set_xticks([0, 1, 2, 3, 4])
cb = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.015)
cb.set_label("Vp (m/s)")
fig.savefig("fig_models.png")
plt.close(fig)

# Convergence
fig, axes = plt.subplots(1, 3, figsize=(10, 3.2), gridspec_kw=dict(wspace=0.32))
for run, name, color, ls in runs:
    o = results[run]
    axes[0].plot(np.arange(len(o["error"])), o["error"], ls, color=color, lw=1.5, label=name)
    axes[2].semilogy(np.arange(len(o["misfit"])), o["misfit"], ls, color=color, lw=1.5)
    if "refl_error" in o:
        axes[1].plot(np.arange(1, len(o["refl_error"]) + 1), o["refl_error"].mean(axis=1), ls, color=color, lw=1.5)
for ax, label, s in zip(axes, ["Vp RMS error (m/s)", "Reflector |error| (m)", "Normalized misfit"], "abc"):
    ax.set_xlabel("Iteration")
    ax.set_ylabel(label)
    ax.text(0.0, 1.02, f"({s})", transform=ax.transAxes, ha="left", va="bottom", fontweight="bold")
axes[1].set_ylim(0, None)
fig.legend(*axes[0].get_legend_handles_labels(), loc="lower center", bbox_to_anchor=(0.5, 0.95), ncol=2, frameon=False)
fig.savefig("fig_convergence.png")
plt.close(fig)

# Traveltime residuals of the source at x = 1.9 km, after the last iteration
xr = np.arange(0.0, 4001.0, 20.0)/1e3
sid, sx = 10, 1.9
obs = np.fromfile(f"data/shot_{sid}_traveltime_p.bin", dtype=np.float32).reshape(nrefl + 1, -1).T
fig, axes = plt.subplots(2, 2, figsize=(10, 5.4), sharex=True, gridspec_kw=dict(hspace=0.28, wspace=0.2))
for run, name, color, ls in runs:
    syn = np.fromfile(f"test_{run}/iteration_{results[run]['last']}/synthetic/shot_{sid}_traveltime_p.bin",
                      dtype=np.float32).reshape(nrefl + 1, -1).T
    for e, ax in enumerate(axes.flat):
        ax.plot(xr, (syn[:, e] - obs[:, e])*1e3, ls, color=color, lw=1.2, label=name)
for e, ax in enumerate(axes.flat):
    ax.axhline(0, color="black", lw=0.6)
    for xb in (sx - 1.0, sx + 1.0):
        ax.axvline(xb, color="black", lw=0.6, ls=":")
    ax.set_title("First arrivals" if e == 0 else f"Reflector {e}", fontsize=10)
    ax.text(0.02, 0.95, f"({'abcd'[e]})", transform=ax.transAxes, ha="left", va="top", fontweight="bold")
for ax in axes[:, 0]:
    ax.set_ylabel("Residual (ms)")
for ax in axes[-1, :]:
    ax.set_xlabel("Receiver x (km)")
fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="lower center", bbox_to_anchor=(0.5, 0.95), ncol=2,
           frameon=False)
fig.savefig("fig_traveltimes.png")
plt.close(fig)

print(f"Vp RMS error above the deepest reflector, from {results['trtt']['error'][0]:.1f} m/s at the start:")
for run, name, *_ in runs:
    print(f"  {name:42s} {results[run]['error'][-1]:5.1f} m/s, normalized misfit {results[run]['misfit'][-1]:.2e}")
for run in ("trtt_all_offsets", "trtt"):
    e = results[run]["refl_error"][-1]
    print(f"{run}: reflector mean |error| " + ", ".join(f"{v:.1f}" for v in e) + " m")
print("Figures: fig_models.png, fig_convergence.png, fig_traveltimes.png")
