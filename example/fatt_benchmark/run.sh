
bindir=$HOME/src/latte/bin

#make clean
#make

# generate model and geometry
#./exec

export OMP_NUM_THREADS=10

## forward modeling
#mpirun -np 1 $bindir/x_eikonal2 param_eikonal.rb

# fatt
mpirun -np 1 $bindir/x_fatt2 param_fatt.rb


