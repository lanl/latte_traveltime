#!/bin/bash

set -e

bindir=../../bin

# Create the models and the geometry; the program needs FLIT and RGM
make
./exec

# One MPI rank per source
export OMP_NUM_THREADS=1

# Compute the traveltimes in the true model
mpirun -np 16 $bindir/x_eikonal3 param_eikonal.rb

# FATT from the initial model
mpirun -np 16 $bindir/x_fatt3 param_fatt.rb

# Plot the models and the data misfit
python3 plot.py
