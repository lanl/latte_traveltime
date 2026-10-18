# Level ranges of the 3D fast sweeps

This test checks a fix of the 3D fast-sweeping solvers. A sweep visits the grid level by level, where a level holds the grid points at the same distance, in grid points, from the starting corner of the sweep. The range of j in a level used nx + ny where it needs nx + nz. On a grid with ny < nz, every sweep therefore skipped some grid points. The other sweeps still reached these points, so the solver still converged, but it needed more iterations.

LATTE stops sweeping when the mean absolute change between two iterations falls below `sweep_stop_threshold` (default `1.0e-4`), or after `sweep_niter_max` iterations (default 10). On the grid of this test, these default settings stopped the solver too early before the fix: its traveltimes differed by up to 256 µs from the fully converged ones. After the fix, they differ by at most 0.5 µs. The fix applies to the forward and the adjoint solvers, and it changes no result when ny ≥ nz.

## Method

The model has 41 × 17 × 31 grid points (x, y, z) at 20 m. The velocity is v = 2000 + 1.5 z m/s with a slow Gaussian anomaly. Three surface sources each record 189 surface receivers. `x_eikonal3` computes the traveltime fields twice: with the default stopping settings, and converged fully (`sweep_stop_threshold = 1.0e-12`, `sweep_niter_max = 200`). The test passes when the two fields agree within 10 µs everywhere. When `x_eikonal3_cuda` exists, its GPU sweeps visit every point of a level, and its fields must equal those of `x_eikonal3` bit for bit. The test then also fails when that run finds no GPU.

## Running the test

Run `./run.sh`, or `bindir=<bin> np=<ranks> ./run.sh` for another build. The test takes a few seconds and needs `python3` with `numpy`. It prints the largest difference and PASS or FAIL, and the exit status is 1 on failure.

## Results

The largest difference is 0.5 µs, and the fields of `x_eikonal3_cuda` equal those of `x_eikonal3` bit for bit.
