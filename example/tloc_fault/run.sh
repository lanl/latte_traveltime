#!/bin/bash

set -e

bindir=../../bin

# Create the models and the geometry; the program needs FLIT and RGM
make clean
make
./exec1

# forward modeling
mpirun -np 40 $bindir/x_eikonal3 param_eikonal.rb

# add noise
./exec2

# tloc without regularization
mpirun -np 40 $bindir/x_tloc3 param.rb

# tloc with smoothness regularization
mpirun -np 40 $bindir/x_tloc3 param_smooth.rb

# tloc with regularization
mpirun -np 40 $bindir/x_tloc3 param_reg.rb
