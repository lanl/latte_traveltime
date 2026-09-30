
nx = 101
nz = 101
dx = 20
dz = 20

ns = 40
file_geometry = geometry/geometry.txt

model_update = vp
file_vp = model/vp_init.bin

process_grad = smooth
grad_smoothx = 40
grad_smoothz = 40

dir_record = data_noisy

niter_max = 60
step_max_vp = 50
min_vp = 2000
max_vp = 4000

model_regularization_method = tgpv
reg_tv_lambda2 = 0
reg_scale_vp = 0.5

dir_working = test_tv
