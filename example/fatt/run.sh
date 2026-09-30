#!/bin/bash

set -e

bindir=../../bin

# Create the models and the geometry; the program needs FLIT and RGM
make clean
make
./exec

export OMP_NUM_THREADS=2

# forward modeling
mpirun -np 10 $bindir/x_eikonal2 param_eikonal.rb

# ad fatt
mpirun -np 10 $bindir/x_fatt2 param_fatt_ad.rb

# dd fatt
mpirun -np 10 $bindir/x_fatt2 param_fatt_dd.rb
