
nx = 401
nz = 201
dx = 10
dz = 10

ns = 20
file_geometry = geometry/geometry.txt

# elastic reflection - s incident
which_medium = elastic-iso
incident_wave = s
model_name = vp, vs, refl
file_vp = model/vp.bin
file_vs = model/vs.bin
file_refl = model/refl.bin
dir_synthetic = data_elastic_refl_s
snaps = 0, 1, 1
dir_snapshot = snapshot_elastic_refl_s
