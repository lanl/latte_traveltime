# 3D FATT with topography: Mount Rainier

This example runs 3D traveltime modeling and first-arrival traveltime tomography (FATT) on a GPU, in a model with the topography of Mount Rainier, Washington.

## Model and geometry

`dem/rainier.tif` is a digital elevation model (DEM) from the 3D Elevation Program of the U.S. Geological Survey, whose data are in the public domain. `dem_to_topo.py` reprojects the DEM to UTM zone 10N, smooths it with a Gaussian of 120 m, and samples it on the model grid. The model covers 12 km by 12 km around the summit with a grid spacing of 60 m, which gives 201 × 201 grid points along x (east) and y (north). It extends 6 km downward from the summit, at an elevation of 4365 m, with 101 grid points in depth. The elevation of the surface ranges from 1388 to 4365 m.

`create_test.f90` builds the P-wave velocity model with RGM. The model simplifies the geology of Mount Rainier. A volcanic edifice of layered lava flows and breccias, whose velocity increases with depth from 3000 to 4200 m/s, overlies an unconformity at an elevation of about 2 km. Below it, layered and faulted volcanic rocks have velocities from 4200 to 5600 m/s. A granodiorite pluton with a velocity of 5900 m/s lies south of the summit, below an elevation of 1300 m. Hydrothermal alteration lowers the velocity beneath the summit by up to 15%, and weathering lowers it near the surface. Below the surface, the velocity ranges from 2400 to 5900 m/s. The initial model is the true model without the pluton and the alteration, with its slowness smoothed over 1.2 km horizontally and 180 m in depth, in coordinates that follow the surface.

The velocity is zero above the surface. The mask `model/mask.bin` is zero there and one at and below the surface. Thirty-six sources lie on a 6 × 6 grid, 2040 m apart, from 0.9 to 11.1 km in x and y. Every source records 1521 receivers on a 39 × 39 grid, 300 m apart, from 0.3 to 11.7 km in x and y. The sources and the receivers lie at the first grid point at or below the surface.

## Running the example

`run.sh` creates the topography, the models and the geometry, computes the traveltimes in the true model with `param_eikonal.rb`, and runs 10 iterations of FATT from the initial model with `param_fatt.rb`. Both runs stop the fast sweeping when the mean change between two iterations falls below `sweep_stop_threshold = 1.0e-6`. FATT masks the gradient above the surface before and after smoothing it over 600 m horizontally and 300 m vertically (`process_grad = mask, smooth, mask`). It limits the velocity change of each iteration to 150 m/s and keeps the velocity within 2000–6500 m/s.

The eikonal solvers run on the GPU with `x_eikonal3_cuda` and `x_fatt3_cuda`, so LATTE must be built with `use_cuda = on` (see `src/README.md`). `run.sh` uses 4 MPI ranks, which share the GPU on a machine with one GPU, with 2 OpenMP threads each, and it expects the LATTE binaries in `../../bin`. `create_test.f90` needs FLIT and RGM, and `dem_to_topo.py` runs with `python3` and needs `rasterio` and `scipy`. The traveltime computation also saves the traveltime field of every source in `snapshot` (`snaps` and `dir_snapshot` in `param_eikonal.rb`), about 590 MB in all, for display.
