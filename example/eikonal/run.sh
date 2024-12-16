
bindir=$HOME/src/latte/bin

make clean
make

# generate model and geometry
./exec

# forward modeling
mpirun -np 20 $bindir/x_eikonal2 param_acoustic.rb
mpirun -np 20 $bindir/x_eikonal2 param_elastic.rb
mpirun -np 20 $bindir/x_eikonal2 param_acoustic_reflection.rb
mpirun -np 20 $bindir/x_eikonal2 param_elastic_reflection_p.rb
mpirun -np 20 $bindir/x_eikonal2 param_elastic_reflection_s.rb

