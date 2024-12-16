
nx = 201
nz = 201
dx = 10
dz = 10

ns = 1
file_geometry = geometry/geometry.txt

model_update = vp
file_vp = model/vp_homo.bin

process_grad = smooth
grad_smoothx = 30
grad_smoothz = 30

niter_max = 1

dir_record = data_low
dir_working = test_low

dir_record = data_high
dir_working = test_high
