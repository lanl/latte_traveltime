# Dirichlet condition against adjoint sources at the receivers

This test compares the adjoint fields of FATT for `adjoint_receiver = nearest` and `adjoint_receiver = linear` with brute-force finite-difference sensitivity kernels, in 2D and 3D. The drivers `kernel2d` and `kernel3d` call `forward_iso` and `adjoint_iso` of LATTE in the same process, so the reference costs one forward solve per perturbation instead of one program run.

## Method

The reference kernel of a receiver is ∂T/∂v_k, the derivative of the receiver traveltime of `forward_iso` with respect to the velocity at grid point k. The driver perturbs the velocity at each grid point by ±1% and takes central differences. The fast sweeping runs until it converges fully. The adjoint kernel is −λ_k/v_k³ times the cell area (2D) or volume (3D), where λ is the adjoint field of a unit residual. For a group of receivers with residuals e_r, the reference gradient of J = ½ Σ e_r² is Σ e_r ∂T_r/∂v_k. The adjoint gradient comes from one adjoint solve with all receivers of the group.

The errors are relative L2 norms. A scaled error first multiplies the adjoint kernel by its least-squares factor. `nearest` needs a factor of about v/dx, which the step-size normalization of FATT absorbs. The cosine is the normalized inner product with the reference. As a check of the reference, the kernels satisfy Σ_k v_k ∂T/∂v_k = −T, because the traveltime scales with 1/v.

## Cases

`2d`: a 121 × 61 grid with 10 m spacing, v = 1800 + 1.2 z m/s with a slow anomaly, and two shots.

- Shot 1 at (100, 0) m. One receiver moves from a grid point across two cells along x, and another along z, in 1 m steps. Two receiver pairs lie in one cell with residuals of +1 ms and −1 ms: one pair rounds to the same grid point, the other to neighboring grid points. A surface line has 44 receivers, 23.7 m apart, with the residuals of a model with two more anomalies.
- Shot 2 at (640, 0) m. A borehole has 60 receivers, 9 m apart, at x = 703 m, with the residuals of the same model.

`3d`: a 41 × 41 × 31 grid with 10 m spacing and one source. One receiver moves along a diagonal across cell faces in x, y and z, and the same two kinds of pairs. The reference covers 11 × 11 × 11 grid points around the receivers.

## Running the test

Build the drivers against the objects of a LATTE build, e.g., `make`, or `make latte2=<objdir of src/2d> latte3=<objdir of src/3d>` for another build. Then run `np=8 ./run.sh`. With 8 MPI ranks, the test takes about 6 minutes, mostly for the 3D reference. The scripts need `python3` with `numpy`, `scipy` and `matplotlib`.

`2d/analyze.py`, `2d/groups.py` and `3d/analyze.py` print the errors. `2d/plot.py` saves `fig_setup.png`, `fig_kernels.png`, `fig_sweep.png`, `fig_pairs.png`, `fig_vsp.png` and `fig_surface.png` in `2d`.

## Results

Errors more than 50 m from the source:

| Case | `linear` | `nearest` |
|---|---|---|
| 2D moving receivers, within 40 m of the receiver | cosine ≥ 0.9999, error 1.8–2.5% | cosine down to 0.74, scaled error up to 27% (x) and 67% (z) in the middle of a cell |
| 2D pair, same grid point | cosine 0.9998, error 2.4% | zero gradient |
| 2D pair, neighboring grid points | cosine 0.9997, error 3.0% | cosine 0.29, scaled error 96% |
| 2D surface line | cosine 0.995, error 10% | cosine 0.977, scaled error 22% |
| 2D borehole | cosine 0.990, error 16% | cosine 0.916, scaled error 40% |
| 3D moving receiver, in the window | cosine 0.99998, error 3.3–3.5% | cosine down to 0.81, scaled error up to 59% |
| 3D pairs | cosine ≥ 0.9999, error 3.5% | zero gradient; cosine 0.85, scaled error 52% |

The amplitude of `linear` does not depend on the receiver position. Its least-squares factor is 0.97 at every position, so its kernel is about 3% larger than the reference near the receiver. The factor of `nearest` changes by up to 11% within a cell in 2D and by 60% in 3D, and by a factor of 2.5 between the borehole receivers. In both treatments, the adjoint kernel is 5–30% too large within 100 m of the source in 2D. The cause is the transport scheme of the adjoint solver, which both treatments share.
