#!/bin/bash

set -e

bindir=../../bin

export OMP_NUM_THREADS=1

# Create the models and the geometry
python3 create_model_and_geometry.py

# Compute the traveltimes in the true model, then add noise to them
mpirun -np 8 $bindir/x_eikonal2 param_eikonal.rb
python3 add_noise.py

# FATT without regularization
mpirun -np 8 $bindir/x_fatt2 param_fatt.rb

# FATT with TGpV regularization
mpirun -np 8 $bindir/x_fatt2 param_fatt_tgpv.rb

# FATT with TV regularization
mpirun -np 8 $bindir/x_fatt2 param_fatt_tv.rb

# FATT with TGpV regularization from iteration 11
mpirun -np 8 $bindir/x_fatt2 param_fatt_tgpv_delayed.rb

# Plot the models and the model errors
python3 plot.py
