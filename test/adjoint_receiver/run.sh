#!/bin/bash

# Compare the adjoint fields of adjoint_receiver = nearest and linear with
# brute-force finite-difference sensitivity kernels, in 2D and in 3D.
# Build the drivers first with make (see the Makefile). np is the number
# of MPI ranks, e.g., np=8 ./run.sh.

set -e

np=${np:-4}
export OMP_NUM_THREADS=1

# ---------------------------------------------------------------- 2D
cd 2d
python3 setup.py
mkdir -p out
dims="121 61 10 10"

# Finite-difference kernels of all receivers, one shot at a time
mpirun -np $np ../kernel2d fd model/v.bin $dims geometry/shot1_all.txt 0.01 out/fd_shot1.bin > out/fd_shot1.log 2>&1
mpirun -np $np ../kernel2d fd model/v.bin $dims geometry/shot2_all.txt 0.01 out/fd_shot2.bin > out/fd_shot2.log 2>&1

# Adjoint fields of each receiver alone
for m in nearest linear; do
    for s in 1 2; do
        mpirun -np $np ../kernel2d each model/v.bin $dims geometry/shot${s}_all.txt $m out/each_shot${s}_$m.bin \
            > out/each_shot${s}_$m.log 2>&1
    done
done

# Residuals of the receiver groups, and their adjoint fields
../kernel2d fwd model/vtrue.bin $dims geometry/shot1_surface.txt out/obs_shot1_surface.trec > out/obs.log 2>&1
../kernel2d fwd model/vtrue.bin $dims geometry/shot2_vsp.txt out/obs_shot2_vsp.trec >> out/obs.log 2>&1
python3 residuals.py
for m in nearest linear; do
    for g in surface pair_same pair_adj; do
        ../kernel2d adj model/v.bin $dims geometry/shot1_$g.txt out/res_$g.bin $m out/adj_${g}_$m.bin \
            > out/adj_${g}_$m.log 2>&1
    done
    ../kernel2d adj model/v.bin $dims geometry/shot2_vsp.txt out/res_vsp.bin $m out/adj_vsp_$m.bin \
        > out/adj_vsp_$m.log 2>&1
done

python3 analyze.py
python3 groups.py
python3 plot.py
cd ..

# ---------------------------------------------------------------- 3D
cd 3d
python3 setup.py
mkdir -p out
dims="41 41 31 10"

# Finite-difference kernels in a window of 11 x 11 x 11 grid points around the receivers
mpirun -np $np ../kernel3d fd model/v.bin $dims geometry/all.txt 0.01 24 34 21 31 13 23 out/fd.bin > out/fd.log 2>&1

for m in nearest linear; do
    mpirun -np $np ../kernel3d each model/v.bin $dims geometry/all.txt $m out/each_$m.bin > out/each_$m.log 2>&1
    for p in pair_same pair_adj; do
        ../kernel3d adj model/v.bin $dims geometry/$p.txt geometry/res_$p.bin $m out/adj_${p}_$m.bin \
            > out/adj_${p}_$m.log 2>&1
    done
done

python3 analyze.py
cd ..
