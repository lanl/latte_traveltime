
nx = 501
nz = 301
dx = 10
dz = 10

ns = 1200
file_geometry = geometry/geometry_no_t0.txt

dir_record = data_noisy

step_max_sx = 1:500, 30:100
step_max_sz = 1:500, 30:100

min_sx = 900.0
max_sx = 4100.0
min_sz = 1450
max_sz = 3000

model_update = sx, sz
model_aux = vp
file_vp = model/vp_init.bin
file_sx = model/sx_init.bin
file_sz = model/sz_init.bin

reg_scale_sx = 0:0, 15:0, 25:0.5, 40:0.8
reg_scale_sz = 0:0, 15:0, 25:0.5, 40:0.8
source_regularization_method = ml
reg_ml_src_infer = $HOME/src/latte/ml/main2_infer.py
reg_ml_model_infer = $HOME/src/latte/ml/infer2.model
reg_ml_src_refine = $HOME/src/latte/ml/main2_refine.py
reg_ml_model_refine = $HOME/src/latte/ml/refine2.model
reg_ml_python = python

misfit_type = dd
dir_working = test_loc_dd_acoustic_reg
