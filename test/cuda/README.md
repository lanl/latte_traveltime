# CPU and CUDA executables of the 3D solvers

This test runs `x_eikonal3`, `x_fatt3`, `x_trtt3` and `x_tloc3` and their CUDA versions, `x_eikonal3_cuda`, `x_fatt3_cuda`, `x_trtt3_cuda` and `x_tloc3_cuda`, on the same parameter files. It times every run and compares the output files of each pair byte by byte. The CUDA kernels repeat the iterations of the CPU solvers, so the results must be the same to the last bit.

## Cases

The model has 101 × 81 × 41 grid points at 12 m, with a velocity gradient, a Gaussian anomaly of +200 m/s and one sinusoidal reflector at about 330 m depth. Four surface sources each record 13 × 11 surface receivers.

- eikonal: the first-arrival and reflection traveltimes in the true model.
- fatt: two iterations from the gradient model, with the first arrivals of the true model.
- trtt: two iterations of `vp` and `refl` from the gradient model, with first arrivals and reflections, and the reflector depths from offsets up to 600 m (`offset_max_refl_depth`).
- tloc: two iterations of `sx`, `sy`, `sz` and `st0` for 8 sources at 300 m depth, from starting positions up to 30 m off, with 7 × 6 surface receivers in the true model.

`x_eikonal3` computes the observed data first.

## Running the test

Build LATTE with `use_cuda = on` (see the README of the repository), and run `./run.sh` on a machine with a GPU, or `bindir=<bin> np=<ranks> ./run.sh` for another build. `python3 compare.py --programs fatt trtt` runs some of the programs only. The test needs `python3` with `numpy`. It prints the times, the speedup and the result of each pair, writes them to `results.txt`, and returns 1 when a pair differs or a CUDA run finds no GPU. The output of every run stays in `<program>_cpu` and `<program>_cuda`, and the logs in `logs`.

## Results

With one GPU and 4 MPI ranks of one thread each, all pairs give identical output files, and the whole test takes about 7 minutes. The times depend on the machine:

| Program | CPU (s) | CUDA (s) | Speedup | Files compared |
|---|---|---|---|---|
| eikonal | 7.4 | 4.2 | 1.8 | 4 |
| fatt | 17.5 | 5.8 | 3.0 | 25 |
| trtt | 145.1 | 62.1 | 2.3 | 32 |
| tloc | 167.1 | 32.0 | 5.2 | 240 |

The four ranks share the GPU, and the grid is small. The speedup grows with the grid: one shot on a 201 × 201 × 161 grid takes 52 s with `x_eikonal3` and 2.4 s with `x_eikonal3_cuda`.
