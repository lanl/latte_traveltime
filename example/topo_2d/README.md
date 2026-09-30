# 2D FATT with topography

This example tests 2D traveltime modeling and first-arrival traveltime tomography (FATT) in a model with a topographic top surface. The velocity is zero above the topography, and FATT masks the gradient there, so the region above the topography stays unchanged. The example also tests TGpV model regularization with this mask.

## Model and geometry

The model is 3 km wide and 1 km deep, with a grid spacing of 10 m (101 × 301 grid points along z and x). `create_test.f90` creates a random topography with 300 m of relief, and a layered and faulted random model with RGM, whose P-wave velocity ranges from 2000 to 4000 m/s. It then increases the velocity by 2000 m/s per kilometer of depth below the topography, so that the first arrivals reach deeper. Below the topography, the velocity ranges from about 2300 to 5700 m/s. The initial model is the true model with a smoothed slowness. Both models have zero velocity above the topography. The mask `model/mask.bin` is zero there and one below.

Thirty sources lie on the topographic surface, 100 m apart. Every source records the receivers at all 301 surface grid points.

## Running the example

`run.sh` compiles and runs `create_test.f90`, computes the traveltimes in the true model, runs 50 iterations of FATT from the initial model without and with regularization, and plots the results. It uses 15 MPI ranks, two shots per rank, and takes a few minutes. The Fortran program needs FLIT and RGM, the Python scripts need `numpy` and `matplotlib`, and `run.sh` expects the LATTE binaries in `../../bin`.

`param_modeling.rb` also saves the full traveltime field of each shot to `snapshot/`. `param_fatt.rb` masks the gradient before and after smoothing it (`process_grad = mask, smooth, mask`): the first mask keeps the region above the topography out of the smoothing, and the second removes what the smoothing spreads into that region.

`plot.py` saves two figures to `figures/`, as PDF and PNG files. `models` shows the true model with the traveltime field of shot 15 and the source positions, the initial model and the inverted model. `misfit` shows the data misfit in each iteration, normalized by its initial value. `plot.py` also prints the root-mean-square (RMS) model error below the topography and the RMS traveltime residual of the initial and inverted models. The model error stays large, because the true model has fine layers and faults that FATT cannot resolve.

## Regularization

`param_fatt_tgpv_0.1.rb`, `param_fatt_tgpv_0.2.rb` and `param_fatt_tgpv_0.5.rb` add TGpV model regularization to `param_fatt.rb`, with strengths (`reg_scale_vp`) of 0.1, 0.2 and 0.5 from iteration 1. `param_fatt_tgpv_delayed.rb` uses the strength 0.2 from iteration 11 (`reg_scale_vp = 1~10:0, 11:0.2`). The denoising also changes the zero velocity above the topography, so the regularization term is nonzero there. LATTE applies the gradient mask again after adding the term, so this region stays unchanged.

`plot_regularization.py` first checks the regularization against the run without regularization:

- LATTE denoises the updated model only when the next iteration has a positive strength, that is, in iterations 1 to 50, or 10 to 50 for the delayed run.
- The iterations before the first regularization term equal those without regularization, bit for bit.
- In the first regularized iteration, the gradient differs from the one without regularization by coef × (m − μ) below the topography. Here m is the model, μ is the denoised model of the previous iteration, and coef = strength × RMS(gradient)/RMS(m − μ).
- The velocity above the topography stays zero.

It then prints the model errors and the traveltime residuals, and saves two figures to `figures/`: `regularization_models` shows the true, initial and inverted models, and `regularization_curves` shows the data misfit and the model error in each iteration. With 15 MPI ranks and one OpenMP thread per rank, we obtain:

| Model | RMS error (m/s) | RMS error, top 50 m (m/s) | RMS residual (ms) |
|---|---|---|---|
| Initial model | 169.2 | 159.3 | 3.92 |
| No regularization | 150.1 | 62.4 | 0.33 |
| TGpV, 0.1 | 150.9 | 74.8 | 0.50 |
| TGpV, 0.2 | 146.2 | 68.1 | 0.39 |
| TGpV, 0.5 | 146.7 | 80.4 | 0.62 |
| TGpV, 0.2 from iteration 11 | 147.8 | 64.2 | 0.35 |

The model errors cover the region below the topography, and the second column its top 50 m. The traveltimes contain no noise, so the regularization raises the residual. The strengths 0.2 and 0.5 lower the model error slightly, but every regularized run has a larger error in the top 50 m, where the denoising meets the zero velocity above the topography. Other numbers of MPI ranks or OpenMP threads change the results slightly.
