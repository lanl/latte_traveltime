
# Introduction

The directory contains scripts for building and running several validation examples. Please refer [our LATTE paper](https://doi.org/10.1093/gji/ggaf079) for details. 

- `eikonal`: Showcasing 2D eikonal equation solving in acoustic and elastic media, for first-arrival and reflection settings. The results are shown in Figures A1-A4 of [the LATTE paper](https://doi.org/10.1093/gji/ggaf079). 
- `eikonal_interpolation`: Comparing linear interpolation used in `LATTE` with conventional nearest-grid interpolation; note that in the code, the nearest-grid interpolation has been disabled, and by default `LATTE` uses linear interpolation. The results are shown in Figures 3 and 4 of [the LATTE paper](https://doi.org/10.1093/gji/ggaf079). 
- `fatt`: An example for validating 2D FATT functionality of `LATTE`. The results are shown in Figures 6-9 of [the LATTE paper](https://doi.org/10.1093/gji/ggaf079). 
- `fatt_benchmark`: A benchmark for validating that `LATTE`'s 2D FATT can generate medium parameter gradients with correct signs. The results are shown in Figure 5 of [the LATTE paper](https://doi.org/10.1093/gji/ggaf079). 
- `fatt_3d`: An example for validating 3D FATT functionality of `LATTE`. 
- `tloc`: An example for validating 2D joint FATT and source location functionality of `LATTE`. The results are shown in Figures 10-21 of [the LATTE paper](https://doi.org/10.1093/gji/ggaf079). 
- `tloc_fault`: An example for validating 3D joint FATT and source location functionality of `LATTE`, as well as ML-enhanced source location associated with faults. The results are shown in Figures 26-36 of [the LATTE paper](https://doi.org/10.1093/gji/ggaf079). 
- `tloc_fracture`: An example for validating 2D source location functionlity of `LATTE`, as well as ML-enhanced source location associated with faults. The results are shown in Figures 22-25 of [the LATTE paper](https://doi.org/10.1093/gji/ggaf079). 
- `yilmaz_near_surface`: An example field dataset (the data only contains picked first-arrival traveltime, not original seismic waveform data and a 1D gradient model is created as the initial Vp model) associated with [this paper (Yilmaz et al., 2022)](
https://doi.org/10.1190/tle41010040.1).

# Misc

To reproduce these results, you need to install several dependencies: 

- [FLIT](https://github.com/lanl/flit) to create the models and source-receiver geometry files with the Fortran codes in the subdirectories. You can also use your own tools to generate these files. 
- [RGM](https://github.com/lanl/rgm) to generate the random geological models used in the examples. 
- [pymplot](https://github.com/lanl/pymplot) to plot results with the Ruby scripts in the subdirectories. 
- The plotting scripts in some of the examples need the three Python scripts in [python](https://github.com/lanl/latte_traveltime/tree/main/misc/python). Please make links to these files in the subfolders in order to plot the results. 