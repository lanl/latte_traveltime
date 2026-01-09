

nx = 97
nz = 21
dx = 1
dz = 1

which_medium = acoustic-iso

ns = 49
file_geometry = ./geometry/geometry.txt

model_update = vp
file_vp = model/vp_init.bin

process_grad = smooth
grad_smoothx = 3
grad_smoothz = 2

min_vp = 350
max_vp = 5500
niter_max = 100

dir_record = pick

dir_working = test_fatt
