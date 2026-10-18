# 2D TRTT: a velocity anomaly above a reflector

This example inverts first-arrival and reflection traveltimes for the P-wave velocity and the depth of a reflector with `x_trtt2`. It compares the result with first-arrival tomography and with an inversion that knows the true reflector.

## Model and geometry

`create_model_and_geometry.py` builds a grid of 301 × 151 points at 10 m. The starting model is v = 1500 + 0.8 z m/s, where z is the depth in meters. The true model adds a Gaussian anomaly of +300 m/s, centered at x = 1.5 km and z = 0.4 km, with standard deviations of 300 m horizontally and 100 m vertically. One reflector follows z = 900 + 100 sin(2πx/3000) m. The reflector model stores it as the label 1 at the nearest grid point of each column plus its offset from that point in grid cells (see `refl` in [parameters.md](../../doc/parameters.md)). Fifteen sources, 200 m apart, and 151 receivers, 20 m apart, lie on the surface.

## Running the example

`run.sh` creates the models and the geometry, and computes the first-arrival and reflection traveltimes in the true model with `x_eikonal2` and `param_eikonal.rb`. It then runs three inversions of 20 iterations with `x_trtt2`. `param_first_arrivals.rb` inverts the first arrivals only, because it gives the reflections zero weight. `param_trtt.rb` inverts the velocity and the reflector from no reflector: in every iteration, TRTT images the reflector, keeps the new image only when it does not raise the data misfit, and refines its depth by least squares. `param_true_reflectors.rb` inverts the velocity with the true reflector fixed, as the best case. Finally, `plot.py` prints the errors and saves `fig_models.png`, `fig_convergence.png` and `fig_traveltimes.png`.

The example needs `python3` with `numpy` and `matplotlib`, and it expects the LATTE executables in `../../bin`. With 4 MPI ranks, it takes about 3 minutes.

## Results

The RMS velocity errors after 20 iterations, from 43.2 m/s in the starting model, are:

| Inversion | Vp RMS error (m/s) |
|---|---|
| First arrivals only | 27.2 |
| TRTT | 21.8 |
| True reflector fixed | 13.6 |

The TRTT reflector lies within 9.1 m of the true one on average, between x = 0.15 and 2.85 km. TRTT fits the reflection traveltimes as well as the inversion with the true reflector, whereas the first-arrival model misses them by up to 26 ms (`fig_traveltimes.png`). The remaining velocity error comes from the trade-off between the velocity above the reflector and its depth. Within the accuracy of the data, a small change of the reflector compensates for part of the velocity error.
