# 3D TRTT: a velocity anomaly above a curved reflector

This example inverts first-arrival and reflection traveltimes for the P-wave velocity and the depth of a reflector with `x_trtt3`. It is the 3D counterpart of [trtt_2d](../trtt_2d), and it compares the result with first-arrival tomography and with an inversion that knows the true reflector.

## Model and geometry

`create_model_and_geometry.py` builds a grid of 61 × 41 × 41 points (x, y, z) at 10 m. The starting model is v = 1500 + z m/s, where z is the depth in meters. The true model adds a Gaussian anomaly of +200 m/s, centered at x = 300 m, y = 200 m and z = 120 m, with standard deviations of 100 m horizontally and 50 m vertically. One reflector follows z = 260 + 30 sin(2πx/600) + 15 cos(2πy/400) m. Twenty-four sources, 100 m apart, and 31 × 21 receivers, 20 m apart, lie on the surface.

## Running the example

`run.sh` creates the models and the geometry, and computes the first-arrival and reflection traveltimes in the true model with `x_eikonal3` and `param_eikonal.rb`. It then runs three inversions of 6 iterations with `x_trtt3`. `param_first_arrivals.rb` inverts the first arrivals only, because it gives the reflections zero weight. `param_trtt.rb` inverts the velocity and the reflector from no reflector. `param_true_reflectors.rb` inverts the velocity with the true reflector fixed. Finally, `plot.py` prints the errors and saves `fig_models.png`, `fig_reflector.png`, `fig_convergence.png` and `fig_traveltimes.png`.

With a build that has `use_cuda = on`, `suffix=_cuda ./run.sh` runs `x_eikonal3_cuda` and `x_trtt3_cuda`, which give the same results. The example needs `python3` with `numpy` and `matplotlib`, and it expects the LATTE executables in `../../bin`. With 4 MPI ranks, it takes about 12 minutes with the CUDA executables, and longer with the CPU ones.

## Results

The RMS velocity errors after 6 iterations, from 32.9 m/s in the starting model, are:

| Inversion | Vp RMS error (m/s) |
|---|---|
| First arrivals only | 31.5 |
| TRTT | 33.1 |
| True reflector fixed | 12.3 |

TRTT images the reflector within 4.9 m on average in the interior, 50 m from the sides, and its data misfit falls to 18% of the starting one. Its velocity, however, does not improve in 6 iterations. The anomaly lies close above the reflector, and the refined reflector absorbs most of its traveltime effect. The inversion with the true reflector shows how much the reflections could add.
