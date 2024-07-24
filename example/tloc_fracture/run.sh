bindir=$HOME/src/latte/bin

make clean
make

# generate model and geometry
./exec1

# forward modeling
mpirun -np 40 $bindir/x_eikonal2 param_eikonal.rb

# add noise
./exec2

# tloc without regularization
mpirun -np 40 $bindir/x_tloc2 param_tloc.rb

# tloc with regularization
mpirun -np 40 $bindir/x_tloc2 param_tloc_reg.rb
