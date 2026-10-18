# 2D TRTT: a faulted layered model

This example inverts first-arrival and reflection traveltimes for the P-wave velocity and the depths of three reflectors with `x_trtt2`, in a layered model with four normal faults. The offsets reach ten times the depth of the shallowest reflector. The example shows that TRTT then needs the reflector depths from the near offsets (`offset_max_refl_depth`). It compares TRTT with first-arrival tomography and with an inversion that knows the true reflectors.

## Model and geometry

`create_model_and_geometry.py` builds the true model with the Python interface of [RGM](https://github.com/lanl/rgm) on a grid of 401 × 161 points at 10 m: twelve layers between 1800 and 3600 m/s, cut by four normal faults. The three reflectors are horizons of the relative geologic time at the largest velocity jumps of the central column, between 0.3 and 1.3 km depth, so they step across the faults. The starting model is the linear fit of the laterally averaged true velocity. Twenty sources, 200 m apart, and 201 receivers, 20 m apart, lie on the surface.

## Running the example

`run.sh` creates the models and the geometry, and computes the first-arrival and reflection traveltimes in the true model with `x_eikonal2` and `param_eikonal.rb`. It then runs four inversions of 20 iterations with `x_trtt2`. `param_first_arrivals.rb` inverts the first arrivals only, because it gives the reflections zero weight. `param_trtt_all_offsets.rb` inverts the velocity and the reflectors from no reflectors, with the reflectors from all offsets. `param_trtt.rb` does the same with `offset_max_refl_depth = 1000`: the reflectors come from the offsets up to 1 km, about twice the depth of the shallowest reflector, while the velocity update uses all offsets. `param_true_reflectors.rb` inverts the velocity with the true reflectors fixed. Finally, `plot.py` prints the errors and saves `fig_models.png`, `fig_convergence.png` and `fig_traveltimes.png`.

The example needs `python3` with `numpy`, `matplotlib` and the Python interface of RGM (see `python/README.md` of RGM), and it expects the LATTE executables in `../../bin`. With 4 MPI ranks, it takes about 20 minutes.

## Results

The RMS velocity errors above the deepest reflector, between x = 0.2 and 3.8 km, after 20 iterations, from 82.8 m/s in the starting model, are:

| Inversion | Vp RMS error (m/s) | Mean reflector depth error (m) |
|---|---|---|
| First arrivals only | 72.4 | – |
| TRTT, reflectors from all offsets | 89.4 | 45.8, 31.5, 32.5 |
| TRTT, reflectors from offsets up to 1 km | 57.1 | 8.5, 7.4, 8.7 |
| True reflectors fixed | 52.7 | – |

With all offsets, the far offsets image the reflectors too deep in the starting model, and the velocity update then compensates for the misplaced reflectors: the velocity error grows, and the misfit stalls at 11% of its starting value. With the reflector depths from the offsets up to 1 km, TRTT fits the data as well as the inversion with the true reflectors and comes within 4.4 m/s of its velocity error. The TRTT reflectors smooth the steps at the faults (`fig_models.png`).
