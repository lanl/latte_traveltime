
nx = 141
ny = 121
nz = 101
dx = 10
dy = 10
dz = 10

ns = 1
file_geometry = geometry/geometry3_gradient.txt

model_name = vp
file_vp = model/v3_gradient.bin

sweep_stop_threshold = 1.0e-10

snaps = 0, 1, 1
dir_snapshot = snapshot3_gradient
dir_synthetic = data3_gradient_interp

#dir_synthetic = data3_gradient_nearest
