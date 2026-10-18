"""Compare the adjoint kernels of both receiver modes with the brute-force
finite-difference kernels. The adjoint gradient of a receiver traveltime with
respect to the velocity at grid point k is -lambda_k/v_k^3*dx*dz."""

import json
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lattetools import read_array

meta = json.load(open("meta.json"))
nx, nz, d = meta["nx"], meta["nz"], meta["d"]
v = read_array("model/v.bin", (nz, nx)).astype(float)
x = np.arange(nx)*d
z = np.arange(nz)*d
Z, X = np.meshgrid(z, x, indexing="ij")


def load_stack(path, nr):
    return read_array(path, (nr, nz, nx)).astype(float)


def grad(lam):
    return -lam/v**3*d*d


def metrics(g, k, mask):
    """cosine, best scale alpha (alpha*g ~ k), relative error of g and of alpha*g"""
    gm, km = g[mask], k[mask]
    cos = gm @ km/np.sqrt((gm @ gm)*(km @ km))
    alpha = (gm @ km)/(gm @ gm)
    return cos, alpha, np.linalg.norm(gm - km)/np.linalg.norm(km), np.linalg.norm(alpha*gm - km)/np.linalg.norm(km)


def shot_data(shot):
    nr = len(np.fromfile(f"out/fd_shot{shot}.bin.trec", dtype=np.float32))
    trec = np.fromfile(f"out/fd_shot{shot}.bin.trec", dtype=np.float32).astype(float)
    fd = load_stack(f"out/fd_shot{shot}.bin", nr)
    near = load_stack(f"out/each_shot{shot}_nearest.bin", nr)
    lin = load_stack(f"out/each_shot{shot}_linear.bin", nr)
    return trec, fd, near, lin


if __name__ == "__main__":
    out = {}
    for shot, src, names in [(1, meta["src1"], ["sweep_x", "sweep_z", "pair_same", "pair_adj", "surface"]),
                             (2, meta["src2"], ["vsp"])]:
        trec, fd, near, lin = shot_data(shot)
        far = np.hypot(X - src[0], Z - src[1]) > 50
        rec = [p for nm in names for p in (meta["groups"][nm] if nm != "vsp" else meta["vsp"])]
        for r in range(len(trec)):
            k = fd[r]
            gn, gl = grad(near[r]), grad(lin[r])
            row = {"rec": rec[r], "T": trec[r],
                   # Euler identity: sum_k K_k v_k = -T
                   "euler_fd": -(k*v).sum()/trec[r], "euler_lin": -(gl*v).sum()/trec[r],
                   "euler_near": -(gn*v).sum()/trec[r]}
            near_rec = np.hypot(X - rec[r][0], Z - rec[r][1]) <= 40
            for tag, mask in [("all", far), ("rec", near_rec & far)]:
                row[f"near_{tag}"] = metrics(gn, k, mask)
                row[f"lin_{tag}"] = metrics(gl, k, mask)
            # smoothed by a Gaussian of 2 cells, as gradient smoothing would
            ks, gns, gls = (gaussian_filter(a, 2.0) for a in (k, gn, gl))
            row["near_smooth"] = metrics(gns, ks, far)
            row["lin_smooth"] = metrics(gls, ks, far)
            out[f"shot{shot}_{r}"] = row
    json.dump(out, open("out/metrics.json", "w"), indent=1, default=float)

    # Print a summary per group
    idx = meta["index"]
    for shot, names in [(1, ["sweep_x", "sweep_z", "pair_same", "pair_adj", "surface"]), (2, ["vsp"])]:
        for nm in names:
            rows = [out[f"shot{shot}_{r}"] for r in idx[nm]]
            print(f"\n== {nm} ({len(rows)} receivers)")
            for key in ["all", "rec", "smooth"]:
                for mode in ["near", "lin"]:
                    m = np.array([r[f"{mode}_{key}"] for r in rows])
                    print(f"  {mode:4s} {key:6s} cos {m[:, 0].min():.4f}-{m[:, 0].max():.4f}  "
                          f"alpha {m[:, 1].min():9.4g}-{m[:, 1].max():9.4g}  "
                          f"relerr {m[:, 2].min():.3f}-{m[:, 2].max():.3f}  scaled relerr {m[:, 3].min():.3f}-{m[:, 3].max():.3f}")
            for e in ["euler_fd", "euler_lin", "euler_near"]:
                a = np.array([r[e] for r in rows])
                print(f"  {e:10s} {a.min():9.4g} .. {a.max():9.4g}")
