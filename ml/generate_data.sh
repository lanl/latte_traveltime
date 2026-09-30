
export OMP_NUM_THREADS=48

./exec2 ntrain=2000 nvalid=200
./exec3 ntrain=1000 nvalid=100

# You can also use the new version of RGM for curved faults:

#./exec2_v2 ntrain=2000 nvalid=200
#./exec3_v2 ntrain=1000 nvalid=100

