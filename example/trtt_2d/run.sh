#!/bin/bash

# 2D TRTT example: the velocity and a reflector from first arrivals and
# reflection traveltimes, against first arrivals only and the true reflector
# fixed. bindir and np can be set in the environment.

set -e

bindir=${bindir:-../../bin}
np=${np:-4}
export OMP_NUM_THREADS=1

# The models and the geometry
python3 create_model_and_geometry.py

# First-arrival and reflection traveltimes in the true model
mpirun -np $np $bindir/x_eikonal2 param_eikonal.rb > log_eikonal.txt 2>&1

# The three inversions, 20 iterations each
for run in first_arrivals trtt true_reflectors; do
    echo "x_trtt2 param_$run.rb"
    mpirun -np $np $bindir/x_trtt2 param_$run.rb > log_$run.txt 2>&1
done

# Errors and figures
python3 plot.py
