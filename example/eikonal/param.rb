
nx = 401
nz = 201
dx = 10
dz = 10

ns = 20
file_geometry = geometry/geometry.txt


## acoustic first-arrival
#which_medium = acoustic-iso
#model_name = vp
#file_vp = model/vp.bin
#dir_synthetic = data_acoustic
#snaps = 0, 1, 1
#dir_snapshot = snapshot_acoustic



## elastic first-arrival
#which_medium = elastic-iso
#model_name = vp, vs
#file_vp = model/vp.bin
#file_vs = model/vs.bin
#dir_synthetic = data_elastic
#snaps = 0, 1, 1
#dir_snapshot = snapshot_elastic



## acoustic reflection
#which_medium = acoustic-iso
#model_name = vp, refl
#file_vp = model/vp.bin
#file_refl = model/refl.bin
#dir_synthetic = data_acoustic_refl
#snaps = 0, 1, 1
#dir_snapshot = snapshot_acoustic_refl



## elastic reflection - p incident
#which_medium = elastic-iso
#incident_wave = p
#model_name = vp, vs, refl
#file_vp = model/vp.bin
#file_vs = model/vs.bin
#file_refl = model/refl.bin
#dir_synthetic = data_elastic_refl_p
#snaps = 0, 1, 1
#dir_snapshot = snapshot_elastic_refl_p


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
