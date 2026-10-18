"""Small helpers for the adjoint-receiver tests: reading and writing LATTE
arrays. Arrays use LATTE's Fortran shapes: (nz, nx) in 2D and (nz, ny, nx) in
3D, z fastest on disk."""

from pathlib import Path

import numpy as np


def write_array(w, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.asarray(w, dtype=np.float32).T.tofile(path)


def read_array(path, shape):
    w = np.fromfile(path, dtype=np.float32, count=int(np.prod(shape)))
    return np.ascontiguousarray(w.reshape(shape[::-1]).T)

