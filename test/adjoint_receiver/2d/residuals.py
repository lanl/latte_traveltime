"""Residuals of the receiver groups: T(v) - T(vtrue) for the surface line and
the borehole, +1 ms and -1 ms for the two receiver pairs."""
import json
import numpy as np
meta = json.load(open("meta.json"))
idx = meta["index"]
t1 = np.fromfile("out/fd_shot1.bin.trec", dtype=np.float32)
t2 = np.fromfile("out/fd_shot2.bin.trec", dtype=np.float32)
res = {
    "surface": t1[idx["surface"]] - np.fromfile("out/obs_shot1_surface.trec", dtype=np.float32),
    "vsp": t2 - np.fromfile("out/obs_shot2_vsp.trec", dtype=np.float32),
    "pair_same": np.array([1e-3, -1e-3]),
    "pair_adj": np.array([1e-3, -1e-3]),
}
for k, r in res.items():
    np.asarray(r, dtype=np.float32).tofile(f"out/res_{k}.bin")
    print(f"{k:10s} n={len(r):3d} residual range {r.min()*1e3:8.3f} .. {r.max()*1e3:8.3f} ms")
