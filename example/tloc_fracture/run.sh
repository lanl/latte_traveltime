#!/bin/bash

set -e

bindir=../../bin

# Create the models and the geometry; the program needs FLIT and RGM
make clean
make
./exec1

# forward modeling
mpirun -np 20 $bindir/x_eikonal2 param_eikonal.rb

# add noise
./exec2

# tloc without regularization
mpirun -np 20 $bindir/x_tloc2 param_tloc.rb

# tloc with regularization
mpirun -np 20 $bindir/x_tloc2 param_tloc_reg.rb
