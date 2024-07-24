
bindir=$HOME/src/latte/bin

make clean
make

# generate model and geometry
./exec

# forward modeling
mpirun -np 20 $bindir/x_eikonal2 param.rb

