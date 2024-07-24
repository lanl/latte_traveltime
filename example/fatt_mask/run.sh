
bindir=$HOME/src/latte/bin

#make clean
#make
#./exec

export OMP_NUM_THREADS=1

$HOME/intel/mpi/bin/mpirun -np 40 $bindir/x_eikonal2 param_eikonal.rb
$HOME/intel/mpi/bin/mpirun -np 40 $bindir/x_fatt2 param_fatt.rb
$HOME/intel/mpi/bin/mpirun -np 40 $bindir/x_fatt2 param_fatt_dd.rb


$HOME/intel/mpi/bin/mpirun -np 40 $bindir/x_eikonal2 param_eikonal_refl.rb
$HOME/intel/mpi/bin/mpirun -np 40 $bindir/x_trtt2 param_trtt.rb
