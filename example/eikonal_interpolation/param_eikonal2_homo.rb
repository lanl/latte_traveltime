
nx = 101
nz = 101
dx = 10
dz = 10

ns = 1
file_geometry = geometry/geometry.txt

model_name = vp
file_vp = model/v_homo.bin

sweep_stop_threshold = 1.0e-10

snaps = 0, 1, 1
dir_snapshot = snapshot_homo
dir_synthetic = data_homo_interp

#dir_synthetic = data_homo_nearest
