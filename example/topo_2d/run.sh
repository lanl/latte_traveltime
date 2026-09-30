#!/bin/bash

set -e

bindir=../../bin

# Create the models, the topography mask and the geometry; the program needs FLIT and RGM
make
./exec

# Two shots per MPI rank
export OMP_NUM_THREADS=1

# Compute the traveltimes in the true model
mpirun -np 15 $bindir/x_eikonal2 param_modeling.rb

# FATT from the initial model
mpirun -np 15 $bindir/x_fatt2 param_fatt.rb

# FATT with TGpV regularization of strength 0.1, 0.2 and 0.5, and of strength 0.2 from iteration 11
mpirun -np 15 $bindir/x_fatt2 param_fatt_tgpv_0.1.rb
mpirun -np 15 $bindir/x_fatt2 param_fatt_tgpv_0.2.rb
mpirun -np 15 $bindir/x_fatt2 param_fatt_tgpv_0.5.rb
mpirun -np 15 $bindir/x_fatt2 param_fatt_tgpv_delayed.rb

# Plot the models, a traveltime field and the data misfit
python3 plot.py

# Check the regularization and compare the inversions
python3 plot_regularization.py
