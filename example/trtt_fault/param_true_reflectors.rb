nx = 401
nz = 161
dx = 10
dz = 10

ns = 20
file_geometry = geometry/geometry.txt

# The velocity with the true reflectors fixed
model_update = vp
model_aux = refl
file_vp = model/vp_init.bin
file_refl = model/refl.bin
misfit_weight = 1, 1, 1, 1

dir_record = data
niter_max = 20
misfit_type = ad
nrefl = 3

process_grad = smooth
grad_smoothx = 100
grad_smoothz = 20
step_max_vp = 50
min_vp = 1000
max_vp = 5000

dir_working = test_true_reflectors
