# 3D FATT with a refraction survey

This example tests 3D first-arrival traveltime tomography (FATT) with a refraction survey. RGM creates a layered and faulted true model. FATT starts from a model whose velocity increases linearly with depth, and it recovers the lateral velocity changes across the faults.

## Model and geometry

The model is 2 km long in x and y and 0.6 km deep, with a grid spacing of 10 m (61 × 201 × 201 grid points along z, y and x). RGM creates 10 layers cut by three faults. The P-wave velocity ranges from 1000 to 3000 m/s and increases with depth. The initial model increases linearly with depth, from the mean velocity at the top of the true model (1310 m/s) to the mean velocity at its bottom (2746 m/s).

Sixteen sources lie on the surface on a 4 × 4 grid, 500 m apart. Every source records 1681 receivers on the surface, on a 41 × 41 grid, 50 m apart.

## Inversion

`param_fatt.rb` runs nine iterations. It smooths the gradient over 80 m horizontally and 40 m vertically, and it limits the velocity change of each iteration to 100 m/s. The model error stops decreasing after about nine iterations: FATT recovers the fault blocks, but not the thin layers.

## Running the example

`run.sh` compiles and runs `create_model_and_geometry.f90`, computes the traveltimes in the true model, runs FATT and plots the results. The Fortran program needs FLIT and RGM, the Python script needs `numpy` and `matplotlib`, and `run.sh` expects the LATTE binaries in `../../bin`. `run.sh` uses 16 MPI ranks, one for each source. With them, the example takes about 30 minutes on a workstation.

`plot.py` saves three figures to `figures/`, as PDF and PNG files. `sections` shows vertical sections through the center of the true, initial and inverted models. `depth_slices` shows the true and inverted models at depths of 0.1, 0.25 and 0.4 km. `misfit` shows the data misfit in each iteration, normalized by its initial value.

`plot.py` also prints the root-mean-square (RMS) model error and traveltime residual of the initial and inverted models. With 16 MPI ranks and one OpenMP thread per rank, we obtain:

| Model | RMS model error (m/s) | RMS traveltime residual (ms) |
|---|---|---|
| Initial | 158.2 | 38.02 |
| FATT, 9 iterations | 101.2 | 4.73 |

The random model depends on the RGM version. These results come from RGM commit 5d49b1a (March 2026).
