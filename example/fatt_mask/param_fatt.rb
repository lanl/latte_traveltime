
nx = 801
nz = 201
dx = 10
dz = 10

ns = 160
file_geometry = geometry/geometry.txt

dir_record = data

model_update = vp
file_vp = model/vp_init.bin

process_grad = smooth, mask
grad_smoothx = 40
grad_smoothz = 30
grad_mask = model/mask.bin

step_max_vp = 1:200, 20:100
min_vp = 1300
max_vp = 3800

dir_working = test_fatt
