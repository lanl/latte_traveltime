# LATTE Parameter Reference

LATTE reads all settings from a plain-text parameter file, which is the first command-line argument, e.g.:

```bash
mpirun -np 20 x_eikonal2 param_eikonal.rb
mpirun -np 20 x_fatt2 param_fatt.rb
```

Each program has a 2-D executable (suffix `2`) and a 3-D executable (suffix `3`):

| Program | Executables | Task |
|---------|-------------|------|
| eikonal | `x_eikonal2`, `x_eikonal3` | Forward first-arrival and reflection traveltimes |
| fatt | `x_fatt2`, `x_fatt3` | First-arrival traveltime tomography |
| trtt | `x_trtt2`, `x_trtt3` | Joint transmission-reflection traveltime tomography |
| tloc | `x_tloc2`, `x_tloc3` | Source location, optionally joint with tomography |

---

## Parameter File Syntax

Each line holds one `key = value` pair. Keys are case-insensitive. LATTE ignores unrecognized keys without a warning, so a misspelled key silently keeps its default value. When a key appears more than once, the last occurrence wins. Reading stops at a line that contains only `exit`.

A line whose key matches no parameter acts as a comment, for example a line that starts with `#`. Do not put a comment after a value on the same line, because LATTE reads it as part of the value. List items are separated by commas. A logical parameter is true for `y`, `yes`, `t`, `true`, `.true.` or `1`, written in all-lowercase or all-uppercase letters, and false for any other value.

At start-up, LATTE copies the parameter file to `<dir_synthetic>/parameters.eikonal.<date-time>` (eikonal) or `<dir_working>/parameters.<program>.<date-time>` (fatt, trtt, tloc). It then reads only this copy, so editing the original file during a run has no effect.

Parameters marked **(iter)** accept an iteration schedule: a comma-separated list of `k:v` or `k1~k2:v` entries, where `k` is an iteration number and `v` is the value. For example, `step_max_vp = 1:50, 20:20` gives 50 at iteration 1 and 20 from iteration 20 on. Between two entries, LATTE interpolates real values linearly, and integer, logical and string values take the value of the nearest entry. Before the first entry and after the last one, the value stays constant. A plain value applies to all iterations.

Lengths use the unit of the grid spacing and the coordinates. Velocities and times must use consistent units, for example m, m/s and s. A parameter marked **yes** in the Required column has no default. Every other parameter takes the listed default.

---

## Directories

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `dir_synthetic` | string | eikonal: output directory for traveltime data | `./data_synthetic` | no |
| `dir_snapshot` | string | eikonal: output directory for full traveltime fields (see `snaps`); must differ from `dir_synthetic` | `./snapshot` | no |
| `dir_working` | string | fatt, trtt, tloc: working directory for all inversion output | `./test` | no |
| `dir_record` | string | fatt, trtt, tloc: directory of the observed traveltime data | `./data` | no |

When LATTE exchanges sources and receivers (always in tloc; see `yn_exchange_sr`), it writes the exchanged observed data to `<dir_working>/record_processed` and reads them from there.

---

## Model Grid (original)

The original grid describes the model files on disk. A model file is a raw little-endian float32 array without a header. The z index runs fastest: the array has the Fortran shape (nz, nx) in 2-D and (nz, ny, nx) in 3-D.

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `nx` | integer | Number of grid points in x | — | **yes** |
| `ny` | integer | Number of grid points in y; `1` in 2-D | `1` | no |
| `nz` | integer | Number of grid points in z | — | **yes** |
| `dx` | float | Grid spacing in x | — | **yes** |
| `dy` | float | Grid spacing in y | `1.0` | no |
| `dz` | float | Grid spacing in z | — | **yes** |
| `ox` | float | x coordinate of the first grid point | `0.0` | no |
| `oy` | float | y coordinate of the first grid point | `0.0` | no |
| `oz` | float | z coordinate of the first grid point | `0.0` | no |

---

## Model Grid (target)

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `nnx` | integer | Target number of grid points in x | `nx` | no |
| `nny` | integer | Target number of grid points in y | `ny` | no |
| `nnz` | integer | Target number of grid points in z | `nz` | no |
| `ddx` | float | Target grid spacing in x | `dx` | no |
| `ddy` | float | Target grid spacing in y | `dy` | no |
| `ddz` | float | Target grid spacing in z | `dz` | no |
| `oox` | float | Target origin in x | `ox` | no |
| `ooy` | float | Target origin in y | `oy` | no |
| `ooz` | float | Target origin in z | `oz` | no |

---

## Computational Domain

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `xmin` | float | Lower x bound of the computational grid | `oox` | no |
| `xmax` | float | Upper x bound of the computational grid | `oox + (nnx - 1)*ddx` | no |
| `ymin` | float | Lower y bound | `ooy` | no |
| `ymax` | float | Upper y bound | `ooy + (nny - 1)*ddy` | no |
| `zmin` | float | Lower z bound | `ooz` | no |
| `zmax` | float | Upper z bound | `ooz + (nnz - 1)*ddz` | no |

LATTE builds the computational grid from the target spacing and these bounds: `ox = xmin` and `nx = nint((xmax - xmin)/ddx) + 1`, and likewise in y and z. When this grid differs from the original grid, LATTE interpolates every model file onto it. The interpolation is linear, except for `refl`, which takes the nearest grid point. LATTE drops sources and receivers outside the bounds. Output models, such as gradients and updated models, live on the computational grid.

---

## Medium

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `which_medium` | string | `acoustic-iso`: P traveltimes from `vp`. `elastic-iso`: P traveltimes from `vp` and S traveltimes from `vs` | `acoustic-iso` | no |
| `data_name` | string list | Names of the data components; the step size reads the files `t<name>_all.bin` | `p` for acoustic media; `p, s` for elastic media | no |
| `incident_wave` | string | Incident wave of elastic reflection traveltimes, `p` or `s`; eikonal and trtt use it when a `refl` model is present | `p` | no |
| `min_vpvsratio` | float | Lower Vp/Vs bound, enforced after each update of an elastic inversion | `1.1` | no |
| `max_vpvsratio` | float | Upper Vp/Vs bound | `9.0` | no |
| `vpvsratio_smoothx` | float | Gaussian smoothing length of the Vp/Vs ratio in x, applied when LATTE enforces the bounds; `0` disables it | `0.0` | no |
| `vpvsratio_smoothy` | float | Same in y (3-D) | `0.0` | no |
| `vpvsratio_smoothz` | float | Same in z | `0.0` | no |

The solvers compute P traveltimes, and S traveltimes for `elastic-iso`, and write them to files named with `p` and `s`. So `data_name` must keep its default for now; LATTE reads it for media with other components, such as qP, qSV and qSH in TTI media. The values `acoustic-tti` and `elastic-tti` pass the input check, but no solver implements them.

When LATTE enforces the Vp/Vs bounds, it adjusts Vs when Vs is updated, alone or together with Vp. It adjusts Vp when Vp is updated and Vs is fixed. Grid points with zero velocity keep zero velocity in this step.

---

## Geometry

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `file_geometry` | string | Geometry index file; each line names one shot file in the same directory | — | **yes** |
| `ns` | integer | Number of shots; LATTE reads the first `ns` lines of the index file, so `ns` must not exceed the number of lines | `1` | no |

A shot file holds the following entries. Each entry starts on a new line, and blank lines between entries are allowed:

```
shot id
number of source points
x y z t0            (one line per source point)
number of receivers
x y z weight        (one line per receiver)
```

Each source and receiver line needs all four numbers. With a missing number, LATTE silently takes the first number of the next line. Text after the last required number of a line is ignored. The index file entries must not contain spaces or slashes. In 2-D, LATTE ignores the y coordinates.

A shot may hold several source points. The computed traveltime is the first arrival from any of them, where point l starts at its own time t0_l. The receiver weight multiplies the misfit of that receiver's data, and a zero weight removes the receiver. Shot IDs must be unique, and LATTE stops otherwise. The number of MPI ranks must not exceed the number of shots that remain after selection. When LATTE exchanges sources and receivers (always in tloc), it must not exceed the number of unique receiver positions.

---

## Shot Selection

LATTE applies these filters after it reads the geometry. It drops a shot when no source point or no receiver remains.

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `src_index` | integer list | Shot sequence numbers (line numbers in the index file) as `first, step, last`; give all three values | `1, 1, ns` | no |
| `src_select` | integer list | Keep only these sequence numbers; overrides `src_index` | none | no |
| `src_exclude` | integer list | Drop these sequence numbers | none | no |
| `sid_min` | integer | Minimum shot ID; the default drops negative IDs | `0` | no |
| `sid_max` | integer | Maximum shot ID | largest integer | no |
| `sid_select` | integer list | Keep only these shot IDs; overrides `src_index` and `src_select` | none | no |
| `sid_exclude` | integer list | Drop these shot IDs | none | no |
| `sxmin`, `sxmax` | float | Drop source points with x outside this range | `-∞`, `+∞` | no |
| `symin`, `symax` | float | Drop source points with y outside this range | `-∞`, `+∞` | no |
| `szmin`, `szmax` | float | Drop source points with z outside this range | `-∞`, `+∞` | no |

---

## Receiver Selection

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `rec_index` | integer list | Receiver sequence numbers within each shot as `first, step, last`; give all three values | `1, 1, ∞` | no |
| `rec_exclude` | integer list | Drop the receivers with these sequence numbers within each shot | none | no |
| `rxmin`, `rxmax` | float | Drop receivers with x outside this range | `-∞`, `+∞` | no |
| `rymin`, `rymax` | float | Drop receivers with y outside this range | `-∞`, `+∞` | no |
| `rzmin`, `rzmax` | float | Drop receivers with z outside this range | `-∞`, `+∞` | no |
| `offset_min` | float | Drop receivers with offsets below this value | `0.0` | no |
| `offset_max` | float | Drop receivers with offsets above this value | `+∞` | no |
| `yn_exchange_sr` | logical | Exchange sources and receivers by reciprocity; LATTE then keeps only the first point of each shot's source | `.false.`; always on in tloc | no |

The offset is the 3-D distance between a receiver and the mean position of the shot's source points, so it includes depth differences. A dropped receiver stays in the receiver list with zero weight. eikonal writes a zero traveltime for it.

With `yn_exchange_sr` in eikonal, the output follows the exchanged geometry: one file per unique receiver position, with one value per shot.

---

## Model Files

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `file_<name>` | string | Binary file of model parameter `<name>` | `''` | see below |

| Name | Programs | Content |
|------|----------|---------|
| `vp` | all | P velocity |
| `vs` | all, `elastic-iso` | S velocity |
| `refl` | eikonal, trtt | Reflector labels: value `l` marks the grid points of reflector `l` (1, 2, ...), and `0` marks all other points |
| `sx`, `sy`, `sz` | tloc | Source coordinates, one float32 value per shot (`sy` in 3-D) |
| `st0` | tloc | Source origin times, one float32 value per shot |

LATTE stops when the file of a `model_aux` entry is missing. For a name in `model_name` or `model_update`, a missing file gives a zero model and only a warning. Always provide these files.

Grid points with zero velocity lie outside the medium. LATTE sets the traveltime to zero there, skips source points in such cells, and returns a zero traveltime for receivers in such cells. A zero-velocity region can therefore model the air above topography.

In tloc, the starting source coordinates and origin times of the inverted parameters come from `file_sx`, `file_sy`, `file_sz` and `file_st0`, not from the geometry. Each file holds one value per shot, in the order of the selected shots. Without such a file, the starting value is zero.

---

## Forward Modeling

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `model_name` | string list | eikonal: models to load. LATTE uses `vp`, `vs` and `refl`, and computes reflection traveltimes when `refl` is present | — | **yes** (eikonal) |
| `sweep_niter_max` | integer | Maximum number of fast-sweeping iterations, for both the forward and the adjoint solver | `10` | no |
| `sweep_stop_threshold` | float | Stop sweeping when the mean absolute change between two iterations falls below this value | `1.0e-4` | no |
| `snaps` | float list | eikonal: when the sum of the values is positive, save each shot's full traveltime field to `dir_snapshot`; the values themselves play no other role | `-1.0` | no |
| `verbose` | logical | Print the grid range of each shot | `.false.` | no |

eikonal writes `shot_<id>_traveltime_p.bin`, and for `elastic-iso` also `shot_<id>_traveltime_s.bin`, to `dir_synthetic`. Each file holds float32 traveltimes, one per receiver in shot-file order. With a `refl` model, a file holds `nrefl + 1` columns, with the receiver index running fastest. The first column holds first arrivals, and column `l + 1` holds the reflection traveltimes of reflector `l`, where `nrefl` is the largest label in `refl`. For `elastic-iso` with `incident_wave = p`, the `_p` file holds P first arrivals and PP reflections. The `_s` file holds a zero first column and PS reflections. With `incident_wave = s`, the `_s` file holds S first arrivals and SS reflections, and the `_p` file holds a zero first column and SP reflections. A snapshot file holds one traveltime field per column of the data file. The fields cover the computational grid, or the shot's subvolume when the adaptive model range is on.

The observed data of fatt, trtt and tloc use the same file names and format in `dir_record`. A negative observed traveltime marks a missing pick, and the inversion ignores it.

---

## Adaptive Model Range (per-shot subvolume)

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `yn_adpx` | logical | Restrict each shot's computation in x to the extent of its source points and active receivers, padded by `adp_extrax` | `.false.` | no |
| `yn_adpy` | logical | Same in y | `.false.` | no |
| `yn_adpz` | logical | Same in z | `.false.` | no |
| `adp_extrax` | float | Padding on each side in x | `0.0` | no |
| `adp_extray` | float | Padding on each side in y | `0.0` | no |
| `adp_extraz` | float | Padding on each side in z | `0.0` | no |
| `adp_taperx` | float | Length of the Blackman taper that LATTE applies to each shot's gradient in x, at the sides of its subvolume that lie inside the model, before stacking | `0.0` | no |
| `adp_tapery` | float | Same in y | `0.0` | no |
| `adp_taperz` | float | Same in z | `0.0` | no |

The taper removes the sharp edges where a shot's gradient ends inside the model. LATTE tapers only along axes with the adaptive range on, so the taper has no effect when `yn_adpx`, `yn_adpy` and `yn_adpz` are all off. Sides on the model boundary stay untapered, and LATTE limits each taper to half the subvolume. A taper no longer than `adp_extrax` stays within the padding, so it leaves the gradient between the sources and receivers unchanged. Reflector images in trtt are not tapered.

---

## Inversion Control (fatt, trtt, tloc)

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `model_update` | string list | Parameters to invert: `vp` and `vs`; `refl` (trtt); `sx`, `sy` (3-D), `sz` and `st0` (tloc) | `vp` | no |
| `model_aux` | string list | Fixed models that the forward solver needs, for example `vs` when an elastic inversion updates only `vp` | none | no |
| `niter_max` | integer | Maximum number of iterations | `100` | no |
| `misfit_type` | string | `ad` (`absolute-difference`) or `dd` (`double-difference`) | `ad` | no |
| `misfit_threshold` | float | Give zero weight to data whose absolute residual exceeds this value | `+∞` | no |
| `misfit_weight` | float list | fatt, tloc: the first value weights all first-arrival data, P and S alike, and LATTE ignores other values. trtt: exactly `nrefl + 1` weights, for the first arrivals and for each reflector | `1.0`; trtt: `nrefl + 1` ones | no |
| `yn_continue` | logical | Resume after the last iteration recorded in `file_data_misfit`; without that file, start from `resume_from_iter` | `.false.` | no |
| `resume_from_iter` | integer | Iteration to start from when `yn_continue` is false or `file_data_misfit` does not exist; LATTE reads the `updated_<name>.bin` files of the previous iteration | `1` | no |
| `yn_flat_stop` | logical | Stop when three successive iterations give the same misfit | `.false.` | no |

When three successive iterations give the same misfit and `yn_flat_stop` is false, LATTE halves every `step_max_<name>` and raises the default `jumpout_factor` to 1.05 for the remaining iterations.

In tloc with `misfit_type = dd` and without `st0` in `model_update`, LATTE differences the traveltimes of each source between receivers, which removes the unknown origin times.

---

## Model Bounds and Step Limits

For each parameter `<name>` in `model_update`:

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `min_<name>` | float | Lower bound, applied after each update | `vp`, `vs`, `refl`, `st0`: `0.0`; `sx`: `ox`; `sy`: `oy`; `sz`: `oz` | no |
| `max_<name>` | float | Upper bound, applied after each update | `vp`, `vs`: `1.0e5`; `refl`: `1.0e9`; `sx`: `ox + (nx - 1)*dx`, and likewise for `sy` and `sz`; `st0`: `1.0e6` | no |
| `step_max_<name>` | float (iter) | Largest change of `<name>` in one update | `vp`, `vs`: `100.0`; `sx`: `0.1*(nx - 1)*dx`, and likewise for `sy` and `sz`; `st0`: `0.1`; `refl`: `1.0` | no |

LATTE does not apply the bounds of `vp` and `vs` at points where the updated value is zero. A zero-velocity region stays zero only when the search direction vanishes there. To keep such a region fixed, mask the gradient with `process_grad = mask` and `grad_mask`.

---

## Search Direction

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `search_method` | string | `SD` (steepest descent), `CG` (nonlinear conjugate gradient, Polak-Ribière; it restarts with steepest descent when its direction is not a descent direction) or `L-BFGS` (10 stored updates). Aliases: `sd`, `steepest-descent`, `cg`, `conjugate-gradient`, `l-bfgs`, `l-BFGS` | `CG` | no |
| `search_method_<name>` | string | Per-parameter override of `search_method`; `refl` always uses `SD` | `search_method` | no |

LATTE does not check the value. An unrecognized value gives a zero search direction, so the parameter never changes.

---

## Step Size

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `step_size_method` | string | `linear`, `quadratic` or `line-search` | `linear` | no |
| `jumpout_factor` | float (iter) | `linear` and `line-search`: accept a trial step when its misfit is below `jumpout_factor` times the current misfit | `1.0`; `1.05` after three equal misfits | no |
| `yn_enforce_update` | logical (iter) | `linear` only: accept the first trial step even when the misfit increases | `.false.` | no |

`linear` estimates the step from one trial model under a linearized misfit and then halves it until the misfit criterion holds, with at most five trials. `line-search` combines quadratic interpolation with bisection. With `linear` and `line-search`, the model stays unchanged in an iteration where no trial meets the criterion. `quadratic` fits a parabola to the misfits of three step sizes and applies the fitted step without a misfit check. LATTE does not check `step_size_method`, so an unrecognized value skips every model update without an error.

---

## Preconditioning

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `yn_precond` | logical | Divide each shot's gradient by its illumination, which LATTE computes as the adjoint field of unit residuals | `.true.` | no |
| `precond_eps` | float | Stabilization of the division, as a fraction of the maximum illumination | `1.0e-4` | no |

When `yn_precond` is false, LATTE divides the adjoint field by v³, where v is the velocity. The division gives infinite or NaN values at zero-velocity points, and LATTE sets them to zero before processing each shot's gradient.

---

## TRTT-specific

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `nrefl` | integer | Number of reflectors. The observed data files hold `nrefl + 1` columns, and LATTE ignores `refl` labels above `nrefl` | `1` | no |
| `offset_min_refl` | float | Give zero weight to reflection data with offsets below this value | `0.0` | no |
| `offset_max_refl` | float | Give zero weight to reflection data with offsets above this value | `+∞` | no |
| `reflector_imaging_threshold` | float | Reflector imaging marks a grid point as a reflector point when the normalized mismatch between the source-side and receiver-side traveltimes, raised to the fourth power, is below this value | `1.0e-9` | no |
| `reflector_<n>_taper_length` | float (iter) | 2-D: length at each lateral end that LATTE excludes when it fits reflector `n`; must be below `nx*dx/3` | `10*dx` | no |
| `reflector_<n>_smooth_window` | float (iter) | 2-D: LOWESS window for smoothing reflector `n` | `0.2*(nx*dx - 2*taper)` | no |
| `reflector_<n>_depth_min` | float (iter) | Shallowest depth of reflector `n` | `oz` | no |
| `reflector_<n>_depth_max` | float (iter) | Deepest depth of reflector `n` | `oz + (nz - 1)*dz` | no |

The reflector index `<n>` starts at 1. The `reflector_<n>_*` parameters apply only when `refl` is in `model_update`. In 3-D, the taper and window parameters have x and y versions: `reflector_<n>_taper_length_x`, `reflector_<n>_taper_length_y`, `reflector_<n>_smooth_window_x` and `reflector_<n>_smooth_window_y`.

When the starting `refl` model is all zero, trtt first images the reflectors from the observed reflection traveltimes. LATTE processes reflector images with the lists `process_shot_refl` and `process_refl`. These lists take the steps of gradient processing, with parameter names that use `refl` in place of `grad`, for example `refl_smoothx`.

---

## Regularization

LATTE regularizes a parameter in two parts. After each model update, it denoises the model with the listed methods. In the next iteration, it adds the difference between the model and its denoised version, times a coefficient, to the gradient. Model regularization acts on `vp` and `vs`, and source regularization acts on the tloc source positions. When the gradient processing of `vp` or `vs` includes `mask`, LATTE applies that mask again after adding the regularization term, so masked zones stay unchanged.

### Enabling Regularization

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `model_regularization_method` | string list | Denoising methods for `vp` and `vs`, applied in order: `tikhonov`, `smooth`, `tgpv`, `structure`, `vpvs_similarity` | none | no |
| `source_regularization_method` | string list | Methods for tloc source positions, applied in order: `clustering-fit`, `ml` | none | no |
| `reg_scale_<name>` | float (iter) | Regularization strength of `<name>`. The coefficient of the gradient term is `reg_scale_<name>` times RMS(gradient)/RMS(model - denoised model) | `reg_scale` | no |
| `reg_scale` | float (iter) | Regularization strength of every parameter without its own `reg_scale_<name>` | `0.0` | no |
| `const_reg` | logical (iter) | Use the fixed coefficient `reg_lambda_<name>` as the regularization strength in place of `reg_scale_<name>` | `.false.` | no |
| `reg_lambda_<name>` | float (iter) | Fixed coefficient of `<name>` when `const_reg` is true | `reg_lambda` | no |
| `reg_lambda` | float (iter) | Fixed coefficient of every parameter without its own `reg_lambda_<name>` when `const_reg` is true | `0.0` | no |
| `rankx`, `ranky`, `rankz` | integer | MPI subdomains of the `tgpv` and `structure` filters (`ranky` in 3-D) | 2-D: floor(sqrt(N)); 3-D: floor(N^(1/3)), where N is the number of MPI ranks | no |

At the end of an iteration, LATTE denoises the updated model when an updated parameter has a positive strength in the next iteration: `vp` or `vs` for the model methods, and `sx`, `sy` or `sz` for the source methods. The next iteration then adds the gradient term with this denoised model. The first iteration of a run, or of a resumed run, has no denoised model yet, so its gradient has no regularization term. A parameter with zero strength is not regularized. LATTE does not regularize `st0`. It ignores unrecognized method names without a warning.

### Tikhonov

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `reg_tikhonov_lambda` | float (iter) | Weight of the gradient-norm (Tikhonov) denoising | `10.0` | no |

### Smooth

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `reg_smoothx` | float (iter) | Gaussian smoothing length in x for all parameters; a negative value selects `reg_smoothx_<name>` | `-1.0` | no |
| `reg_smoothx_<name>` | float (iter) | Gaussian smoothing length in x for `<name>` | `dx` | no |
| `reg_smoothy`, `reg_smoothy_<name>` | float (iter) | Same in y (3-D) | `-1.0`, `dy` | no |
| `reg_smoothz`, `reg_smoothz_<name>` | float (iter) | Same in z | `-1.0`, `dz` | no |

For `vp` and `vs`, LATTE smooths the slowness 1/v.

### TGpV

TGpV denotes total generalized p-variation. It has a first-order term with weight α₀ and a second-order term with weight α₁, and `reg_tv_norm` sets the power p of both. LATTE has no separate `tv` method: for TV, list `tgpv` and set `reg_tv_lambda2 = 0`, which removes the second-order term. In the defaults below, m is the model to regularize.

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `reg_tv_mu_<name>` | float (iter) | Regularization weight μ for `<name>` | `100/max(abs(m))` | no |
| `reg_tv_lambda1` | float (iter) | Weight α₀ of the first-order term; in 2-D, `-1` selects `reg_tv_lambda1_<name>` | 2-D: `-1.0`; 3-D: `1.0` | no |
| `reg_tv_lambda1_<name>` | float (iter) | 2-D only: α₀ for `<name>` | `1.0` | no |
| `reg_tv_lambda2` | float (iter) | Weight α₁ of the second-order term; `0` gives TV; in 2-D, `-1` selects `reg_tv_lambda2_<name>` | 2-D: `-1.0`; 3-D: `1.0` | no |
| `reg_tv_lambda2_<name>` | float (iter) | 2-D only: α₁ for `<name>` | `1.0` | no |
| `reg_tv_norm` | float (iter) | Power p of the TGpV terms, for example `0.5` or `1.0` | `0.5` | no |
| `reg_tv_niter` | integer (iter) | Number of TGpV iterations | `50` | no |

### Structure

The method `structure` applies anisotropic diffusion filtering (ANDF). Its lengths are in grid points.

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `reg_andf_alpha` | float (iter) | ANDF weight α | `0.001` | no |
| `reg_andf_beta` | float (iter) | ANDF weight β | `1.0` | no |
| `reg_andf_gamma` | float (iter) | 3-D only: ANDF weight γ | `1.0` | no |
| `reg_andf_smoothx` | float (iter) | ANDF smoothing length in x | `2.0` | no |
| `reg_andf_smoothy` | float (iter) | 3-D only: ANDF smoothing length in y | `2.0` | no |
| `reg_andf_smoothz` | float (iter) | ANDF smoothing length in z | `8.0` | no |
| `reg_andf_t` | integer (iter) | Number of diffusion steps | 2-D: `10`; 3-D: `5` | no |
| `reg_andf_sigma` | float (iter) | Diffusion length σ | `10.0` | no |
| `reg_andf_powerm` | float (iter) | Power-law exponent m | `4.0` | no |
| `reg_andf_aux` | string (iter) | Auxiliary model that guides the structure | `''` | no |
| `reg_andf_coh` | string (iter) | Coherence model | `''` | no |

### Vp/Vs Similarity

The method `vpvs_similarity` acts on `vs` only and needs `vp` in `model_update` or `model_aux`. It clips, median-filters and smooths the Vp/Vs ratio, and then sets Vs = Vp/ratio.

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `reg_similarity_smoothx` | float (iter) | Gaussian smoothing length of the ratio in x | `5*dx` | no |
| `reg_similarity_smoothy` | float (iter) | 3-D only: smoothing length in y | `5*dy` | no |
| `reg_similarity_smoothz` | float (iter) | Smoothing length in z | `5*dz` | no |
| `reg_similarity_vpvs_ratio_min` | float (iter) | Lower bound of the ratio | `min_vpvsratio` | no |
| `reg_similarity_vpvs_ratio_max` | float (iter) | Upper bound of the ratio | `max_vpvsratio` | no |

### Clustering-fit Source Regularization (tloc)

The method `clustering-fit` groups the sources with HDBSCAN. It then moves the sources of each cluster onto a fitted curve (2-D) or surface (3-D). In the defaults below, N_s is the number of sources.

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `clustering_fit_min_sample` | integer (iter) | HDBSCAN minimum number of samples | `nint(0.25*N_s)` | no |
| `clustering_fit_min_cluster_size` | integer (iter) | HDBSCAN minimum cluster size | 2-D: `clustering_fit_min_sample`; 3-D: `nint(0.25*N_s)` | no |
| `clustering_fit_method` | string (iter) | Fitting method | `polynomial` | no |
| `clustering_fit_smooth` | float (iter) | Smoothness of the fit | `0.5` | no |
| `clustering_fit_order` | integer (iter) | Polynomial order | 2-D: `1`; 3-D: `2` | no |

### ML Source Regularization (tloc)

The method `ml` rasterizes the current source positions and predicts a fault probability map with a trained network. It then moves each source to the nearest grid point whose probability is at least 0.5. It changes the source coordinates only, not `st0`.

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `reg_ml_python` | string | Python executable | — | **yes** for `ml` |
| `reg_ml_src` | string | Inference script, such as `ml/main2.py` (2-D) or `ml/main3.py` (3-D) | — | **yes** for `ml` |
| `reg_ml_model_infer` | string | Trained inference model | — | **yes** for `ml` |
| `reg_ml_model_refine` | string | Trained refinement model | — | **yes** for `ml` |
| `reg_ml_niter_refine` | integer | Number of refinement iterations | `3` | no |
| `reg_ml_xyz_weight` | float list | Distance weights: `wx, wz` in 2-D, and `wx, wy, wz` in 3-D | all `1.0` | no |
| `reg_ml_max_dist` | float (iter) | Move a source only when its weighted distance to the nearest fault point, sqrt(wx·Δx² + wz·Δz²) (with wy·Δy² added in 3-D), is at most this value; Δx, Δy and Δz are the coordinate differences | `+∞` | no |

---

## Gradient Processing

LATTE processes gradients in two stages. `process_shot_grad` lists the steps for each shot's gradient before stacking, and `process_grad` lists the steps for the stacked gradient. The steps run in the listed order. By default, the same lists apply to every updated model parameter, and `uniform_processing` switches to separate lists for each parameter. LATTE does not process the gradients of `sx`, `sy`, `sz` and `st0`.

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `process_shot_grad` | string list | Steps for each shot's gradient | none | no |
| `process_grad` | string list | Steps for the stacked gradient | none | no |
| `uniform_processing` | logical | Use `process_shot_grad` and `process_grad` for all parameters. When false, each parameter `<name>` uses its own lists, `process_shot_grad_<name>` and `process_grad_<name>`, and its step parameters carry the same infix, for example `shot_grad_vp_smoothx` and `grad_vp_smoothx` | `.true.` | no |
| `<name>_update_iter` | integer list | Update parameter `<name>` only from iteration `first` to `last`; a single value means `first` to `niter_max`. Outside this range, LATTE zeroes the gradient | `1, niter_max` | no |

LATTE ignores unrecognized step names without a warning. In 3-D, each parameter that has x and z versions also has a y version, for example `shot_grad_smoothy`.

### Per-Shot Gradient Processing

Before the listed steps, LATTE sets NaN and infinite values of each shot's gradient to zero. The per-shot parameters do not accept iteration schedules.

| Step | Parameters (defaults) | Description |
|------|-----------------------|-------------|
| `smooth` | `shot_grad_smoothx`, `shot_grad_smoothz` (`3*dx`, `3*dz`) | Gaussian smoothing |
| `maxbal` | — | Divide by the maximum absolute value |
| `rmsbal` | — | Divide by the RMS value |
| `movingbal` | `shot_grad_movingbalx`, `shot_grad_movingbalz` (2-D: `3*dx`, `3*dz`; 3-D: `6*dx`, `6*dy`, `6*dz`) | Moving-window amplitude balancing |
| `medianfilt` | `shot_grad_medianfiltx`, `shot_grad_medianfiltz` (`dx`, `dz`) | Median filtering |
| `dipfilt` | `shot_grad_dipfiltzx` (`-100, 0, 100`), `shot_grad_dipfiltzx_amps` (`0, 0, 0`); 3-D adds the pairs `shot_grad_dipfiltzy` and `shot_grad_dipfiltyx` | Dip filtering; no effect while all amplitudes are zero |
| `removenan` | — | Replace NaN and infinite values with zero |
| `laplacefilt` | — | Laplacian filtering |
| `wavenumfilt` | `shot_grad_wavenumx`, `shot_grad_wavenumx_amps`, `shot_grad_wavenumz`, `shot_grad_wavenumz_amps` (`-1`, which disables the axis) | Wavenumber-domain filtering along each axis |
| `taper` | `shot_grad_taperx`, `shot_grad_taperz`: one length, or `begin, end` (`0, 0`) | Blackman taper at all edges of the shot's subvolume |
| `mask` | `shot_grad_mask` (one file for all shots), or `dir_shot_grad_mask` (a directory with one file `shot_<id>_mask.bin` per shot, which takes precedence); without either, `mask` has no effect | Multiply by a model-size mask, cropped to the shot's range |
| `conemute` | 2-D: `shot_grad_conemutex`, `shot_grad_conemutez` (`-1`, off), `shot_grad_conemutepower` (`2.0`), `shot_grad_conemutetaper` (`10*dx`) | Zero the gradient above the source. Below the source, keep a cone whose half-width grows to `shot_grad_conemutex` at `shot_grad_conemutez` below the source, and taper outside it |

In 3-D, `conemute` takes the radius `shot_grad_conemuter` in place of `shot_grad_conemutex`. It also takes `shot_grad_conemuteorigin`, set to `source` (default) or `receiver`, and its default taper is `10*mean(dx, dy, dz)`.

### Global Gradient Processing

All global step parameters accept iteration schedules, except `grad_taperx`, `grad_tapery`, `grad_taperz` and `grad_andf_rankx`, `grad_andf_ranky`, `grad_andf_rankz`.

| Step | Parameters (defaults) | Description |
|------|-----------------------|-------------|
| `smooth` | `grad_smoothx`, `grad_smoothz` (`3*dx`, `3*dz`) | Gaussian smoothing |
| `maxbal` | — | Divide by the maximum absolute value |
| `rmsbal` | — | Divide by the RMS value |
| `rmsbalx` | `grad_rmsbalx` (`dx`) | Divide the gradient at each x position by its L2 norm over a window of width `grad_rmsbalx` centered there |
| `rmsbalz` | `grad_rmsbalz` (`dz`) | Same along z |
| `rmsbaly`, `rmsbalxy` | `grad_rmsbaly` (`dy`), and `grad_rmsbalx` for `rmsbalxy` | 3-D only: same along y, or over windows in the x-y plane |
| `movingbal` | `grad_movingbalx`, `grad_movingbalz` (2-D: `6*dx`, `6*dz`; 3-D: `3*dx`, `3*dy`, `3*dz`) | Moving-window amplitude balancing |
| `medianfilt` | `grad_medianfiltx`, `grad_medianfiltz` (`dx`, `dz`) | Median filtering |
| `taper` | `grad_taperx`, `grad_taperz`: one length, or `begin, end` (`0, 0`) | Blackman taper at the grid edges |
| `mask` | `grad_mask` (file); without it, `mask` has no effect | Multiply by a model-size mask |
| `scale` | `grad_scale` (`1.0`) | Multiply by a constant |
| `signed_power` | `grad_signed_power` (`0.5`) | Replace each value g by sign(g)·abs(g)^p, where p is `grad_signed_power` |
| `andf` | `grad_andf_*`, see below | Anisotropic diffusion filtering |

The parameters of the `andf` step use physical lengths, unlike `reg_andf_*`:

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `grad_andf_smoothx` | float | ANDF smoothing length in x | `2*dx` |
| `grad_andf_smoothy` | float | 3-D only: ANDF smoothing length in y | `2*dy` |
| `grad_andf_smoothz` | float | ANDF smoothing length in z | `8*dz` |
| `grad_andf_powerm` | float | Power-law exponent | `1.0` |
| `grad_andf_t` | integer | Number of diffusion steps | `5` |
| `grad_andf_sigma` | float | Diffusion length | `6*max(dx, dz)`; 3-D: `6*max(dx, dy, dz)` |
| `grad_andf_alpha` | float | ANDF weight α | `1.0e-3` |
| `grad_andf_beta` | float | ANDF weight β | `1.0` |
| `grad_andf_gamma` | float | 3-D only: ANDF weight γ | `1.0` |
| `grad_andf_aux` | string | Auxiliary model file | `''` |
| `grad_andf_coh` | string | Coherence model file | `''` |
| `grad_andf_rankx`, `grad_andf_ranky`, `grad_andf_rankz` | integer | 3-D only: MPI subdomains of the filter | `1` |

---

## Output Files

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `file_data_misfit` | string | ASCII misfit history; each line holds the iteration number, the misfit and the misfit normalized by the initial misfit | `<dir_working>/data_misfit.txt` | no |
| `file_shot_misfit` | string | float32 array of per-shot misfits with shape (N_s, k + 1) for iterations 0 to k, shot index fastest | `<dir_working>/shot_misfit.bin` | no |

For each iteration k, `<dir_working>/iteration_<k>/model/` holds `grad_<name>.bin`, `srch_<name>.bin` and `updated_<name>.bin`, plus `reg_<name>.bin` when the denoising runs. `<dir_working>/iteration_<k>/synthetic/` holds the synthetic data of that iteration. `iteration_0` holds the starting models and their synthetic data.
