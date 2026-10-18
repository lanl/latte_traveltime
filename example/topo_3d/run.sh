#!/bin/bash

set -e

bindir=../../bin

# The topography from the DEM; the script needs rasterio and scipy
python3 dem_to_topo.py

# The models and the geometry; the program needs FLIT and RGM
make
./exec

# The eikonal solvers run on the GPU with the CUDA executables, so LATTE must be
# built with use_cuda = on (see src/README.md); on a machine with one GPU, the
# four MPI ranks share it
export OMP_NUM_THREADS=2

# Compute the traveltimes in the true model
mpirun -np 4 $bindir/x_eikonal3_cuda param_eikonal.rb

# FATT from the initial model
mpirun -np 4 $bindir/x_fatt3_cuda param_fatt.rb
