
nx = 301
nz = 101

dx = 10
dz = 10

ns = 30
file_geometry = ./geometry/geometry.txt


which_medium = acoustic-iso
model_update = vp

file_vp = model/vp_init.bin

dir_record = data

process_grad = mask, smooth, mask
grad_mask = model/mask.bin
grad_smoothx = 30
grad_smoothz = 20

niter_max = 50

dir_working = test_tgpv_0.2

sweep_stop_threshold = 1.0e-6

verbose = y

model_regularization_method = tgpv
reg_scale_vp = 0.2
