"""Group gradients: the true gradient of J = 1/2 sum_r e_r^2, sum_r e_r K_r from the
finite-difference kernels, against the adjoint gradients of both modes run with
all receivers of the group."""
import json
import numpy as np
from scipy.ndimage import gaussian_filter
from analyze import meta, read_array, grad, metrics, load_stack, nx, nz, X, Z

idx = meta["index"]
res = {}
for g, shot, src in [("pair_same", 1, meta["src1"]), ("pair_adj", 1, meta["src1"]),
                     ("surface", 1, meta["src1"]), ("vsp", 2, meta["src2"])]:
    nr = len(np.fromfile(f"out/fd_shot{shot}.bin.trec", dtype=np.float32))
    fd = load_stack(f"out/fd_shot{shot}.bin", nr)
    e = np.fromfile(f"out/res_{g}.bin", dtype=np.float32).astype(float)
    gfd = np.tensordot(e, fd[idx[g]], axes=1)
    far = np.hypot(X - src[0], Z - src[1]) > 50
    print(f"\n== {g}: |true gradient| max {np.abs(gfd[far]).max():.3e}")
    for m in ["nearest", "linear"]:
        ga = grad(read_array(f"out/adj_{g}_{m}.bin", (nz, nx)).astype(float))
        if np.abs(ga[far]).max() == 0:
            print(f"  {m:8s} zero gradient")
            res[f"{g}_{m}"] = dict(cos=np.nan, alpha=np.nan, relerr=1.0, srelerr=np.nan, cos_s=np.nan, srelerr_s=np.nan)
            continue
        c, a, r0, r1 = metrics(ga, gfd, far)
        cs, as_, _, r1s = metrics(gaussian_filter(ga, 2.0), gaussian_filter(gfd, 2.0), far)
        print(f"  {m:8s} cos {c:7.4f} alpha {a:9.4g} relerr {r0:6.3f} scaled relerr {r1:6.3f} | smoothed: cos {cs:7.4f} scaled relerr {r1s:6.3f}")
        res[f"{g}_{m}"] = dict(cos=c, alpha=a, relerr=r0, srelerr=r1, cos_s=cs, srelerr_s=r1s)
    np.save(f"out/gfd_{g}.npy", gfd)
json.dump(res, open("out/group_metrics.json", "w"), indent=1)
