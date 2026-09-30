
nx = 301
nz = 101

dx = 10
dz = 10

ns = 30
file_geometry = ./geometry/geometry.txt


which_medium = acoustic-iso
model_name = vp

file_vp = model/vp.bin
dir_synthetic = data

#file_vp = model/vp_init.bin
#dir_synthetic = data_init

sweep_stop_threshold = 1.0e-6

snaps = 0, 1, 1

verbose = y
