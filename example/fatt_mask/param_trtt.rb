
nx = 801
nz = 201
dx = 10
dz = 10

ns = 160
file_geometry = geometry/geometry.txt

dir_record = data_refl

nrefl = 1
model_update = vp
model_aux = refl
file_vp = model/vp_init.bin
file_refl = model/refl.bin

process_grad = smooth, mask
grad_smoothx = 40
grad_smoothz = 30
grad_mask = model/mask.bin

step_max_vp = 1:200, 20:100
min_vp = 1300
max_vp = 3800

dir_working = test_trtt

niter_max = 25
misfit_weight = 1, 5
#yn_continue = y
