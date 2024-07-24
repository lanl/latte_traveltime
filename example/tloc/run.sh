
bindir=$HOME/src/latte/bin

make clean
make

# generate model and geometry
./exec

export OMP_NUM_THREADS=1

# forward modeling
mpirun -np 40 $bindir/x_eikonal2 param_eikonal.rb

# ad tloc
mpirun -np 40 $bindir/x_tloc2 param_loc_elastic.rb

# dd tloc
mpirun -np 40 $bindir/x_tloc2 param_loc_dd_elastic.rb

# tomo dd tloc
mpirun -np 40 $bindir/x_tloc2 param_tomo_dd_elastic.rb

# tomo dd tloc with model parameter regularization
mpirun -np 40 $bindir/x_tloc2 param_tomo_dd_elastic_reg.rb

