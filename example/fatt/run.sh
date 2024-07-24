
bindir=$HOME/src/latte/bin

#make clean
#make
#./exec

export OMP_NUM_THREADS=4

mpirun -np 12 $bindir/x_eikonal2 param_eikonal.rb

mpirun -np 12 $bindir/x_fatt2 param_fatt_ad.rb

mpirun -np 12 $bindir/x_fatt2 param_fatt_dd.rb


