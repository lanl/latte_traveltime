bindir=$HOME/src/latte/bin

mpirun -np 40 $bindir/x_eikonal3 param_eikonal.rb
mpirun -np 40 $bindir/x_tloc3 param_tloc.rb
mpirun -np 40 $bindir/x_tloc3 param_tloc_reg.rb
