
bindir=$HOME/src/latte/bin

make clean
make

# generate model and geometry
./exec

export OMP_NUM_THREADS=4

# forward modeling
mpirun -np 10 $bindir/x_eikonal2 param_eikonal.rb

# ad fatt
mpirun -np 10 $bindir/x_fatt2 param_fatt_ad.rb

# dd fatt
mpirun -np 10 $bindir/x_fatt2 param_fatt_dd.rb


