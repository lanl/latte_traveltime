
nx = 501
nz = 301
dx = 10
dz = 10

ns = 300
file_geometry = geometry/geometry_no_t0.txt

which_medium = elastic-iso

dir_record = data

model_update = sx, sz, vp, vs
file_vp = model/vp_init.bin
file_vs = model/vs_init.bin
file_sx = model/sx_init.bin
file_sz = model/sz_init_deep.bin

min_vp = 1600
max_vp = 2600
min_vs = 900
max_vs = 1400
min_sx = 900
max_sx = 4100
min_sz = 50
max_sz = 2950

step_max_sx = 1:500, 20:100
step_max_sz = 1:500, 20:100
step_max_vp = 1:50, 20:20
step_max_vs = 1:50, 20:20

process_grad = smooth
grad_smoothx = 30
grad_smoothz = 30

misfit_type = dd
dir_working = test_tomo_dd_elastic
