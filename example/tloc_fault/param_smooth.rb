
nx = 301
ny = 201
nz = 101
dx = 10
dy = 10
dz = 10

ns = 1200
file_geometry = geometry/geometry_no_t0.txt

which_medium = elastic-iso

step_max_sx = 400
step_max_sy = 400
step_max_sz = 300
step_max_vp = 200
step_max_vs = 200

min_sx = 0.0
max_sx = 3000.0
min_sy = 0.0
max_sy = 2000.0
min_sz = 0.0
max_sz = 1000.0
min_vp = 800
max_vp = 3400
min_vs = 300
max_vs = 2400

process_grad = smooth
grad_smoothx = 30
grad_smoothy = 30
grad_smoothz = 30

model_update = sx, sy, sz, vp, vs
file_vp = model/vp_init.bin
file_vs = model/vs_init.bin
file_sx = model/sx_init.bin
file_sy = model/sy_init.bin
file_sz = model/sz_init.bin

model_regularization_method = smooth
reg_scale_vp = 1:0, 30:0.5
reg_scale_vs = 1:0, 30:0.5
reg_smoothx = 40
reg_smoothy = 40
reg_smoothz = 40

misfit_type = dd

niter_max = 50

dir_record = data_noisy
dir_working = test_tloc_smooth
