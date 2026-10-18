nx = 201
ny = 201
nz = 101
dx = 60
dy = 60
dz = 60

ns = 36
file_geometry = geometry/geometry.txt

which_medium = acoustic-iso
model_update = vp
file_vp = model/vp_init.bin

dir_record = data

process_grad = mask, smooth, mask
grad_mask = model/mask.bin
grad_smoothx = 600
grad_smoothy = 600
grad_smoothz = 300

niter_max = 10
step_max_vp = 150
min_vp = 2000
max_vp = 6500

sweep_stop_threshold = 1.0e-6

dir_working = test

