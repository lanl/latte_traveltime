nx = 61
ny = 41
nz = 41
dx = 10
dy = 10
dz = 10

ns = 24
file_geometry = geometry/geometry.txt

# TRTT: the velocity and the reflectors, from no reflectors
model_update = vp, refl
file_vp = model/vp_init.bin
file_refl = model/refl_zero.bin

dir_record = data
niter_max = 6
misfit_type = ad
nrefl = 1

process_grad = smooth
grad_smoothx = 30
grad_smoothy = 30
grad_smoothz = 20
step_max_vp = 50

dir_working = test_trtt
