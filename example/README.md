
# Introduction

The directory contains scripts for building and running several validation examples. Please refer [our LATTE paper](https://doi.org/10.1093/gji/ggaf079) for details. 

- `eikonal`: Showcasing 2D eikonal equation solving in acoustic and elastic media, for first-arrival and reflection settings. The results are shown in Figures A1-A4 of [the LATTE paper](https://doi.org/10.1093/gji/ggaf079). 
- `eikonal_interpolation`: Comparing linear interpolation used in `LATTE` with conventional nearest-grid interpolation; note that in the code, the nearest-grid interpolation has been disabled, and by default `LATTE` uses linear interpolation. The results are shown in Figures 3 and 4 of [the LATTE paper](https://doi.org/10.1093/gji/ggaf079). 
- `eikonal_multi_source`: Showcasing 2D and 3D eikonal equation solving for a source made of many points, each with its own start time: two point sources, and a rupturing fault represented by point sources along a line (2D) or over a plane (3D). The Python scripts `multi_source_2d.py` and `multi_source_3d.py` build the models and geometries, run `LATTE`, and plot the results; they need `numpy` and `matplotlib`. 
- `fatt`: An example for validating 2D FATT functionality of `LATTE`. The results are shown in Figures 6-9 of [the LATTE paper](https://doi.org/10.1093/gji/ggaf079). 
- `fatt_benchmark`: A benchmark for validating that `LATTE`'s 2D FATT can generate medium parameter gradients with correct signs. The results are shown in Figure 5 of [the LATTE paper](https://doi.org/10.1093/gji/ggaf079). 
- `fatt_3d`: An example for validating 3D FATT functionality of `LATTE`, with a refraction survey over a layered and faulted random model on a 61 × 201 × 201 grid. The Fortran code creates the models and the geometry with [RGM](https://github.com/lanl/rgm), and the Python script `plot.py` plots the results; it needs `numpy` and `matplotlib`. 
- `fatt_regularization`: A checkerboard example for validating the model regularization of 2D FATT in `LATTE`. It inverts noisy traveltimes without regularization, with TGpV or TV regularization, and with TGpV regularization from iteration 11. The Python scripts create the models and the geometry, add noise to the traveltimes, and plot the results; they need `numpy` and `matplotlib`. 
- `tloc`: An example for validating 2D joint FATT and source location functionality of `LATTE`. The results are shown in Figures 10-21 of [the LATTE paper](https://doi.org/10.1093/gji/ggaf079). 
- `tloc_fault`: An example for validating 3D joint FATT and source location functionality of `LATTE`, as well as ML-enhanced source location associated with faults. The results are shown in Figures 26-36 of [the LATTE paper](https://doi.org/10.1093/gji/ggaf079). 
- `tloc_fracture`: An example for validating 2D source location functionlity of `LATTE`, as well as ML-enhanced source location associated with faults. The results are shown in Figures 22-25 of [the LATTE paper](https://doi.org/10.1093/gji/ggaf079). 
- `topo_2d`: An example for validating 2D eikonal equation solving and FATT of `LATTE` for a model with topographic top surface; this is done by masking the region above the topography. It also tests TGpV model regularization with this mask. The Fortran code creates the models and the geometry with [RGM](https://github.com/lanl/rgm), and the Python scripts `plot.py` and `plot_regularization.py` plot the results and check the regularization; they need `numpy` and `matplotlib`. 
- `yilmaz_near_surface`: An example field dataset (the data only contains picked first-arrival traveltime, not original seismic waveform data and a 1D gradient model is created as the initial Vp model) associated with [this paper (Yilmaz et al., 2022)](
https://doi.org/10.1190/tle41010040.1).

# Misc

To reproduce these results, you need to install several dependencies: 

- [FLIT](https://github.com/lanl/flit) to create the models and source-receiver geometry files with the Fortran codes in the subdirectories. You can also use your own tools to generate these files. 
- [RGM](https://github.com/lanl/rgm) to generate the random geological models used in the examples. `fatt`, `tloc_fracture` and `tloc_fault` use the RGM types `rgm2` and `rgm3`, which RGM 2.0 keeps for backward compatibility; the other examples use `rgm2_curved` and `rgm3_curved`. 
- [pymplot](https://github.com/lanl/pymplot) to plot results with the Ruby scripts in the subdirectories. 
- The plotting scripts in some of the examples need the three Python scripts in [python](https://github.com/lanl/latte_traveltime/tree/main/misc/python). Please make links to these files in the subfolders in order to plot the results. 