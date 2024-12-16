
nx = 401
nz = 201
dx = 10
dz = 10

ns = 20
file_geometry = geometry/geometry.txt

# acoustic first-arrival
which_medium = acoustic-iso
model_name = vp
file_vp = model/vp.bin
dir_synthetic = data_acoustic
snaps = 0, 1, 1
dir_snapshot = snapshot_acoustic
