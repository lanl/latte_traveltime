
nx = 121
nz = 61
dx = 10
dz = 10

ns = 1

model_update = vp
file_vp = ../../model/v.bin

# One iteration: iteration_1/model/grad_vp.bin holds the gradient. Without
# preconditioning and processing, it is the adjoint field divided by -v^3.
niter_max = 1
yn_precond = n
step_max_vp = 1
min_vp = 1000
max_vp = 3000

sweep_stop_threshold = 1.0e-12
sweep_niter_max = 100

