#
# Add Gaussian noise with a standard deviation of 3 ms to the traveltimes
#
import numpy as np
from pathlib import Path

ns = 40

rng = np.random.default_rng(2026)
Path('data_noisy').mkdir(exist_ok=True)
for i in range(1, ns + 1):
    t = np.fromfile(f'data/shot_{i}_traveltime_p.bin', np.float32).astype(np.float64)
    # Keep the traveltimes near the sources positive
    t = np.maximum(t + rng.normal(0.0, 3.0e-3, t.size), 1.0e-4)
    t.astype(np.float32).tofile(f'data_noisy/shot_{i}_traveltime_p.bin')
