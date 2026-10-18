"""3D: compare the adjoint kernels of both receiver modes with the
finite-difference kernels inside the FD window."""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lattetools import read_array

meta = json.load(open("meta.json"))
nx, ny, nz, d = meta["nx"], meta["ny"], meta["nz"], meta["d"]
v = read_array("model/v.bin", (nz, ny, nx)).astype(float)
w = meta["window"]
win = np.zeros((nz, ny, nx), bool)
win[w[4] - 1:w[5], w[2] - 1:w[3], w[0] - 1:w[1]] = True
nr = sum(len(g) for g in meta["groups"].values())
trec = np.fromfile("out/fd.bin.trec", dtype=np.float32).astype(float)


def grad(lam):
    return -lam/v**3*d**3


def metrics(g, k, mask):
    gm, km = g[mask], k[mask]
    cos = gm @ km/np.sqrt((gm @ gm)*(km @ km))
    alpha = (gm @ km)/(gm @ gm)
    return dict(cos=cos, alpha=alpha, relerr=np.linalg.norm(gm - km)/np.linalg.norm(km),
                srelerr=np.linalg.norm(alpha*gm - km)/np.linalg.norm(km))


fd = read_array("out/fd.bin", (nr, nz, ny, nx)).astype(float)
lam = {m: read_array(f"out/each_{m}.bin", (nr, nz, ny, nx)).astype(float) for m in ["nearest", "linear"]}
out = {"receivers": [], "pairs": {}}
idx = meta["index"]
print("sweep receivers (window metrics; Euler = -sum(g v)/T over the whole grid)")
for r in range(nr):
    row = {"T": trec[r]}
    for m in ["nearest", "linear"]:
        g = grad(lam[m][r])
        row[m] = metrics(g, fd[r], win)
        row[m]["euler"] = -(g*v).sum()/trec[r]
    out["receivers"].append(row)
    if r in idx["sweep"]:
        print(f"  s={2*idx['sweep'].index(r):2d}  nearest cos {row['nearest']['cos']:.4f} srelerr {row['nearest']['srelerr']:.3f} "
              f"alpha {row['nearest']['alpha']:7.2f} | linear cos {row['linear']['cos']:.5f} relerr {row['linear']['relerr']:.4f} "
              f"srelerr {row['linear']['srelerr']:.4f} alpha {row['linear']['alpha']:.4f} euler {row['linear']['euler']:.4f}")
for p in ["pair_same", "pair_adj"]:
    e = np.fromfile(f"geometry/res_{p}.bin", dtype=np.float32).astype(float)
    gfd = np.tensordot(e, fd[idx[p]], axes=1)
    print(f"{p}: |true gradient| max {np.abs(gfd[win]).max():.3e}")
    for m in ["nearest", "linear"]:
        g = grad(read_array(f"out/adj_{p}_{m}.bin", (nz, ny, nx)).astype(float))
        if np.abs(g[win]).max() == 0:
            print(f"  {m:8s} zero gradient")
            out["pairs"][f"{p}_{m}"] = None
            continue
        mm = metrics(g, gfd, win)
        out["pairs"][f"{p}_{m}"] = mm
        print(f"  {m:8s} cos {mm['cos']:.4f} alpha {mm['alpha']:.4g} relerr {mm['relerr']:.3f} srelerr {mm['srelerr']:.3f}")
json.dump(out, open("out/metrics3d.json", "w"), indent=1, default=float)
