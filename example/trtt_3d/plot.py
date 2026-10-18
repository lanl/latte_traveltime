"""Errors and figures of the 3D TRTT example."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({"font.family": "Arial", "font.size": 10, "savefig.dpi": 300, "savefig.bbox": "tight",
                     "savefig.pad_inches": 0.02})
nx, ny, nz, d = 61, 41, 41, 10.0
x, y, z = np.arange(nx)*d, np.arange(ny)*d, np.arange(nz)*d
interior = (slice(5, 36), slice(5, 56))          # columns 50 m from the sides
runs = [("first_arrivals", "First arrivals only", "blue", "-."),
        ("trtt", "TRTT", "red", "-"),
        ("true_reflectors", "True reflector fixed", "black", ":")]


def load(path):
    return np.fromfile(path, dtype=np.float32).reshape(nx, ny, nz).transpose(2, 1, 0)


def iterations(run):
    return sorted(int(p.name.split("_")[1]) for p in Path(f"test_{run}").glob("iteration_*")
                  if (p/"model"/"updated_vp.bin").is_file())


def depth_of(refl, label=1):
    """The depth of a reflector in each (y, x) column of a label-plus-offset model."""
    iz = np.floor(refl + 0.5)
    zz = np.where(iz == label, (np.arange(nz)[:, None, None] + refl - label)*d, np.nan)
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
        out["refl_error"] = [float(np.nanmean(np.abs(depth_of(r) - depth)[interior])) for r in refl]
    results[run] = out

# Velocity change: sections at y = 200 m and slices at z = 120 m, and the reflector
jy, kz = 20, 12
fig, axes = plt.subplots(2, 4, figsize=(10, 4.6), gridspec_kw=dict(hspace=0.32, wspace=0.08, height_ratios=[0.67, 1.0]))
panels = [("True model", vp_true, None)] + [(name, results[run]["vp"], run) for run, name, *_ in runs]
for col, (name, v, run) in enumerate(panels):
    ax = axes[0, col]
    im = ax.imshow((v - vp_init)[:, jy, :], extent=[0, x[-1]/1e3, z[-1]/1e3, 0], cmap="seismic", vmin=-200, vmax=200,
                   aspect="equal", interpolation="none")
    ax.plot(x/1e3, depth[jy]/1e3, "--" if run == "trtt" else "-", color="black", lw=1.0)
    if run == "trtt":
        ax.plot(x/1e3, results[run]["depth"][jy]/1e3, "-", color="black", lw=1.2)
    ax.set_title(name if run is None else f"{name}\n{results[run]['error'][-1]:.1f} m/s", fontsize=9.5)
    ax.set_xlabel("x (km)")
    bx = axes[1, col]
    bx.imshow((v - vp_init)[kz], extent=[0, x[-1]/1e3, y[-1]/1e3, 0], cmap="seismic", vmin=-200, vmax=200,
              aspect="equal", interpolation="none")
    bx.set_xlabel("x (km)")
    for a, s in ((ax, "abcd"[col]), (bx, "efgh"[col])):
        a.text(0.03, 0.95, f"({s})", transform=a.transAxes, ha="left", va="top", fontweight="bold")
        if col > 0:
            a.set_yticklabels([])
axes[0, 0].set_ylabel("z (km)\nat y = 0.2 km")
axes[1, 0].set_ylabel("y (km)\nat z = 0.12 km")
cb = fig.colorbar(im, ax=axes, fraction=0.02, pad=0.015)
cb.set_label("Vp change (m/s)")
fig.savefig("fig_models.png")
plt.close(fig)

# The reflector: true depth and TRTT depth error
fig, axes = plt.subplots(1, 2, figsize=(10, 3.0), gridspec_kw=dict(wspace=0.3))
im0 = axes[0].imshow(depth, extent=[0, x[-1]/1e3, y[-1]/1e3, 0], cmap="viridis", aspect="equal")
im1 = axes[1].imshow(results["trtt"]["depth"] - depth, extent=[0, x[-1]/1e3, y[-1]/1e3, 0], cmap="seismic",
                     vmin=-15, vmax=15, aspect="equal", interpolation="none")
fig.colorbar(im0, ax=axes[0], fraction=0.035, pad=0.02).set_label("True depth (m)")
fig.colorbar(im1, ax=axes[1], fraction=0.035, pad=0.02).set_label("TRTT depth error (m)")
for ax, s in zip(axes, "ab"):
    ax.set_xlabel("x (km)")
    ax.set_ylabel("y (km)")
    ax.text(0.03, 0.95, f"({s})", transform=ax.transAxes, ha="left", va="top", fontweight="bold",
            bbox=dict(facecolor="white", edgecolor="none", pad=0.6))
fig.savefig("fig_reflector.png")
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

# Traveltime residuals along the receiver line y = 200 m, after the last iteration
xr = np.arange(0.0, 601.0, 20.0)/1e3
line = np.arange(31)*21 + 10                      # receivers at y = 200 m
shots = {10: (0.25, 0.15), 19: (0.45, 0.25)}
fig, axes = plt.subplots(2, 2, figsize=(10, 5.0), sharex=True, gridspec_kw=dict(hspace=0.12, wspace=0.25))
for col, (sid, (sx, sy)) in enumerate(shots.items()):
    obs = np.fromfile(f"data/shot_{sid}_traveltime_p.bin", dtype=np.float32).reshape(2, -1).T
    for run, name, color, ls in runs:
        syn = np.fromfile(f"test_{run}/iteration_{results[run]['last']}/synthetic/shot_{sid}_traveltime_p.bin",
                          dtype=np.float32).reshape(2, -1).T
        for e in range(2):
            axes[e, col].plot(xr, (syn[line, e] - obs[line, e])*1e3, ls, color=color, lw=1.2, label=name)
    axes[0, col].set_title(f"Source at (x, y) = ({sx:g}, {sy:g}) km", fontsize=10)
    for e in range(2):
        axes[e, col].axhline(0, color="black", lw=0.6)
        axes[e, col].text(0.03, 0.95, f"({'abcd'[2*e + col]})", transform=axes[e, col].transAxes, ha="left",
                          va="top", fontweight="bold")
    axes[1, col].set_xlabel("Receiver x (km) at y = 0.2 km")
axes[0, 0].set_ylabel("First arrival\nresidual (ms)")
axes[1, 0].set_ylabel("Reflection\nresidual (ms)")
fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="lower center", bbox_to_anchor=(0.5, 0.93), ncol=3,
           frameon=False)
fig.savefig("fig_traveltimes.png")
plt.close(fig)

print(f"Vp RMS error, from {results['trtt']['error'][0]:.1f} m/s at the start:")
for run, name, *_ in runs:
    print(f"  {name:22s} {results[run]['error'][-1]:5.1f} m/s, normalized misfit {results[run]['misfit'][-1]:.2e}")
print(f"TRTT reflector mean |error| in the interior: {results['trtt']['refl_error'][-1]:.1f} m")
print("Figures: fig_models.png, fig_reflector.png, fig_convergence.png, fig_traveltimes.png")
