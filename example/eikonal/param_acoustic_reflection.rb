
nx = 401
nz = 201
dx = 10
dz = 10

ns = 20
file_geometry = geometry/geometry.txt

# acoustic reflection
which_medium = acoustic-iso
model_name = vp, refl
file_vp = model/vp.bin
file_refl = model/refl.bin
dir_synthetic = data_acoustic_refl
snaps = 0, 1, 1
dir_snapshot = snapshot_acoustic_refl
