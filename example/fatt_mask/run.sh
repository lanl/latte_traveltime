
bindir=$HOME/src/latte/bin

make clean
make

# generate model and geometry
./exec

export OMP_NUM_THREADS=1

# fatt
mpirun -np 40 $bindir/x_eikonal2 param_eikonal.rb
mpirun -np 40 $bindir/x_fatt2 param_fatt.rb
mpirun -np 40 $bindir/x_fatt2 param_fatt_dd.rb

# reflection tomography
mpirun -np 40 $bindir/x_eikonal2 param_eikonal_refl.rb
mpirun -np 40 $bindir/x_trtt2 param_trtt.rb
