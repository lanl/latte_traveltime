
nx = 401
nz = 51
dx = 10
dz = 10

ns = 40
file_geometry = geometry/geometry.txt

model_update = vp
file_vp = model/vp_init.bin

process_grad = smooth
grad_smoothx = 30
grad_smoothz = 30

dir_record = data

niter_max = 100
step_max_vp = 50
min_vp = 480
max_vp = 2500

misfit_type = dd
dir_working = test_dd

