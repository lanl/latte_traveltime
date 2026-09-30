# FATT with model regularization

This example tests the model regularization of 2D first-arrival traveltime tomography (FATT) on a checkerboard model. The observed traveltimes contain random noise. Without regularization, the inversion fits the noise, and the model error grows in late iterations. With TGpV regularization, the inverted model stays close to the checkerboard. TV regularization flattens the model into large blocks, the staircase effect of TV, and its final error exceeds that of the inversion without regularization.

## Model and data

The model is 2 km wide and 2 km deep, with a grid spacing of 20 m (101 by 101 grid points). The true model is a checkerboard of 5 by 5 cells, each 400 m wide, with P-wave velocities of 2850 and 3150 m/s. The initial model has a constant velocity of 3000 m/s. Forty sources lie on the four sides of the model, 10 on each side. Two hundred receivers lie on the four sides, 40 m apart. We add Gaussian noise with a standard deviation of 3 ms to the traveltimes.

## Inversions

All four inversions run 60 iterations and differ only in the regularization:

| Parameter file | Regularization | Output directory |
|---|---|---|
| `param_fatt.rb` | None | `test_noreg` |
| `param_fatt_tgpv.rb` | TGpV from iteration 1 | `test_tgpv` |
| `param_fatt_tv.rb` | TV from iteration 1 | `test_tv` |
| `param_fatt_tgpv_delayed.rb` | TGpV from iteration 11 | `test_tgpv_delayed` |

TGpV denotes total generalized _p_-variation. It has first- and second-order terms, and their power _p_ is 0.5 by default. Without the second-order term, TGpV becomes TV (total variation), which `param_fatt_tv.rb` selects with `reg_tv_lambda2 = 0`. With `reg_scale_vp = 0.5`, the regularization term has half the root-mean-square amplitude of the data-misfit gradient.

In `param_fatt_tgpv_delayed.rb`, the schedule `reg_scale_vp = 1~10:0, 11:0.5` sets the strength to zero in iterations 1 to 10. These iterations have no regularization term, so they equal those of `test_noreg`. LATTE also skips the denoising in these iterations, except at the end of iteration 10, where it prepares the regularized model for iteration 11.

## Running the example

`run.sh` creates the models and the geometry, computes the traveltimes, adds the noise, runs the four inversions and plots the results. It needs the LATTE binaries in `../../bin`, and `numpy` and `matplotlib` for the Python scripts. With 8 MPI ranks, it takes about a minute on a workstation.

`plot.py` saves the figure as `figures/checkerboard_regularization.pdf` and `.png`. The figure shows (a) the true model, (b–e) the models after 60 iterations, and (f) the model error in each iteration. The model error is the root-mean-square difference between the inverted and true models. `plot.py` also prints the errors. With 8 MPI ranks and one OpenMP thread per rank, we obtain:

| Inversion | Final error (m/s) | Minimum error (m/s) | Iteration of minimum |
|---|---|---|---|
| No regularization | 79.8 | 73.7 | 27 |
| TGpV | 72.1 | 71.4 | 49 |
| TV | 98.3 | 91.7 | 28 |
| TGpV from iteration 11 | 73.3 | 71.8 | 37 |

Other numbers of MPI ranks or OpenMP threads change the results slightly.
