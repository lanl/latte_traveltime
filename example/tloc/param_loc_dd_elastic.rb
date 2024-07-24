
nx = 501
nz = 301
dx = 10
dz = 10

ns = 300
file_geometry = geometry/geometry_no_t0.txt

which_medium = elastic-iso

dir_record = data

step_max_sx = 1:500, 10:100
step_max_sz = 1:500, 10:100

min_sx = 900
max_sx = 4100
min_sz = 50
max_sz = 2950

model_update = sx, sz
model_aux = vp, vs
file_vp = model/vp.bin
file_vs = model/vs.bin
file_sx = model/sx_init.bin
file_sz = model/sz_init.bin

misfit_type = dd
dir_working = test_loc_dd_elastic
