
nx = 801
nz = 201
dx = 10
dy = 10
dz = 10

ns = 160
file_geometry = geometry/geometry.txt

which_medium = acoustic-iso
nrefl = 1
model_name = vp, refl
file_vp = model/vp.bin
file_refl = model/refl.bin

dir_synthetic = data_refl

snaps = 0,1,1
dir_snapshot = snapshot_refl
