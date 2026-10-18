#!/bin/bash

# 3D TRTT example: the velocity and a reflector from first arrivals and
# reflection traveltimes, against first arrivals only and the true reflector
# fixed. bindir and np can be set in the environment; suffix=_cuda runs the
# CUDA executables of a build with use_cuda = on, which give the same results.

set -e

bindir=${bindir:-../../bin}
np=${np:-4}
suffix=${suffix:-}
export OMP_NUM_THREADS=1

# The models and the geometry
python3 create_model_and_geometry.py

# First-arrival and reflection traveltimes in the true model
mpirun -np $np $bindir/x_eikonal3$suffix param_eikonal.rb > log_eikonal.txt 2>&1

# The three inversions, 6 iterations each
for run in first_arrivals trtt true_reflectors; do
    echo "x_trtt3$suffix param_$run.rb"
    mpirun -np $np $bindir/x_trtt3$suffix param_$run.rb > log_$run.txt 2>&1
done

# Errors and figures
python3 plot.py
