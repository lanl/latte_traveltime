nx = 401
nz = 161
dx = 10
dz = 10

ns = 20
file_geometry = geometry/geometry.txt

# TRTT: the velocity and the reflectors, from no reflectors; the reflectors
# from all offsets
model_update = vp, refl
file_vp = model/vp_init.bin
file_refl = model/refl_zero.bin

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

dir_working = test_trtt_all_offsets
