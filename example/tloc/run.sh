
bindir=$HOME/src/__latte/bin

#make clean
#make
#./exec
#

export OMP_NUM_THREADS=5

$HOME/intel/mpi/bin/mpirun -np 50 $bindir/x_eikonal2 param_eikonal.rb

$HOME/intel/mpi/bin/mpirun -np 50 $bindir/x_tloc2 param_loc_elastic.rb

$HOME/intel/mpi/bin/mpirun -np 50 $bindir/x_tloc2 param_loc_dd_elastic.rb

$HOME/intel/mpi/bin/mpirun -np 50 $bindir/x_tloc2 param_tomo_dd_elastic.rb

$HOME/intel/mpi/bin/mpirun -np 50 $bindir/x_tloc2 param_tomo_dd_elastic_reg.rb

