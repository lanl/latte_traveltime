
nx = 301
ny = 201
nz = 101
dx = 10
dy = 10
dz = 10

ns = 1200
file_geometry = geometry/geometry_no_t0.txt

which_medium = elastic-iso

step_max_sx = 400
step_max_sy = 400
step_max_sz = 300
step_max_vp = 200
step_max_vs = 200

min_sx = 0.0
max_sx = 3000.0
min_sy = 0.0
max_sy = 2000.0
min_sz = 0.0
max_sz = 1000.0
min_vp = 800
max_vp = 3400
min_vs = 300
max_vs = 2400

process_grad = smooth
grad_smoothx = 30
grad_smoothy = 30
grad_smoothz = 30

model_update = sx, sy, sz, vp, vs
file_vp = model/vp_init.bin
file_vs = model/vs_init.bin
file_sx = model/sx_init.bin
file_sy = model/sy_init.bin
file_sz = model/sz_init.bin

model_regularization_method = tgpv, vpvs_similarity
reg_scale_vp = 1:0, 15:0, 30:0.5
reg_scale_vs = 1:0, 15:0, 30:0.5
reg_tv_mu_vp = 0.1
reg_tv_mu_vs = 0.15
reg_tv_lambda1 = 1
reg_tv_lambda2 = 1
reg_similarity_smoothx = 50
reg_similarity_smoothy = 50
reg_similarity_smoothz = 50
reg_similarity_vpvs_ratio_min = 1.3
reg_similarity_vpvs_ratio_max = 2.3

source_regularization_method = ml
reg_scale_sx = 1:0, 15:0, 30:0.5
reg_scale_sy = 1:0, 15:0, 30:0.5
reg_scale_sz = 1:0, 15:0, 30:0.5
reg_ml_src_infer = /users/kaigao/src/latte/ml/main3_infer.py
reg_ml_model_infer = /users/kaigao/src/latte/ml/infer3.model
reg_ml_src_refine = /users/kaigao/src/latte/ml/main3_refine.py
reg_ml_model_refine = /users/kaigao/src/latte/ml/refine3.model
reg_ml_python = /users/kaigao/.conda/envs/my_root/bin/python

misfit_type = dd

niter_max = 50
#yn_continue = y

dir_record = data_noisy
dir_working = test_tloc_reg
