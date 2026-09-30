
nx = 201
ny = 201
nz = 61
dx = 10
dy = 10
dz = 10

ns = 16
file_geometry = geometry/geometry.txt

model_update = vp
file_vp = model/vp_init.bin

process_grad = smooth
grad_smoothx = 30
grad_smoothy = 30
grad_smoothz = 10

dir_record = data

niter_max = 9
step_max_vp = 100
min_vp = 500
max_vp = 4000

dir_working = test
