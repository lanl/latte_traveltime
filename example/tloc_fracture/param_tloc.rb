
nx = 501
nz = 301
dx = 10
dz = 10

ns = 1200
file_geometry = geometry/geometry_no_t0.txt

dir_record = data_noisy

step_max_sx = 1:500, 30:100
step_max_sz = 1:500, 30:100

min_sx = 900.0
max_sx = 4100.0
min_sz = 1450
max_sz = 3000

model_update = sx, sz
model_aux = vp
file_vp = model/vp_init.bin
file_sx = model/sx_init.bin
file_sz = model/sz_init.bin

niter_max = 50
misfit_type = dd
dir_working = test_loc_dd_acoustic
