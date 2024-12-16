
nx = 401
nz = 201
dx = 10
dz = 10

ns = 20
file_geometry = geometry/geometry.txt

# elastic first-arrival
which_medium = elastic-iso
model_name = vp, vs
file_vp = model/vp.bin
file_vs = model/vs.bin
dir_synthetic = data_elastic
snaps = 0, 1, 1
dir_snapshot = snapshot_elastic
