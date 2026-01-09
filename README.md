# Description
**LATTE: Los Alamos TravelTime package based on Eikonal equation**

`LATTE` is an open-source software for 
- **Traveltime computation** based on eikonal equation in 2D/3D acoustic/elastic media
    - First-arrival P, S traveltimes
    - Both first-arrival and reflection traveltimes (e.g., P, S, PP, SS, PS, SP) for an arbitrary number of reflectors. 
- **Adjoint-state first-arrival traveltime tomography** based on first-arrival traveltime in 2D/3D acoustic/elastic media
- **Source location** based on first-arrival traveltime in 2D/3D acoustic/elastic media, including
    - Source location/relocation using P, S, or both P and S traveltimes
    - Joint source location/relocation with first-arrival traveltime tomography using P, S, or both P and S traveltimes

Algorithm features (details are explained in [the LATTE paper](https://academic.oup.com/gji/article/241/2/1275/8046728)):

- Mostly written in modern Fortran. 
- Applies to Cartesian, regularly sampled grid. 
- Forward and adjoint-state equations are solved with the fast-sweeping method. 
- Traveltime computation is based on the factorized eikonal equation for avoiding source singularity and thus improved accuracy. 
- Two-level parallelization, including shot parallelization based on MPI and shared-memory parallelization for each shot based on OpenMP. 
- Various inversion schemes, including 
    - SD: steepest descent
    - NCG: nonlinear conjugate gradient
    - _l_-BFGS: limited-memory Broyden-Fletcher-Goldfarb-Shanno
- Various regularization schemes, including 
    - TV: total variation regularization for medium parameters
    - TGpV: total generalized _p_-variation regularization for medium parameters
    - smooth: Tikhonov-like smoothing regularization for medium parameters
    - ML: machine-learning-based regularization for source location parameters
- User-friendly parameter setting based on [FLIT](https://github.com/lanl/flit) flexible parameter input functionality. 
- User-friendly geometry setting. 

`LATTE` does not yet support:
- Curvilinear or unstructured mesh. 
- Medium anisotropy.
- GPU-based acceleration. 

The work is supported by Los Alamos National Laboratory (LANL) Laboratory Directory Research and Development (LDRD) project 20240322ER. LANL is operated by Triad National Security, LLC, for the National Nuclear Security Administration (NNSA) of the U.S. Department of Energy (DOE) under Contract No. 89233218CNA000001. The research uses high-performance computing resources provided by LANL's Institutional Computing (IC) program. 

The codes are released under LANL open source approval reference O4770.

# Requirement
`LATTE` depends on [FLIT](https://github.com/lanl/flit). Some examples in [example](example) use [RGM](https://github.com/lanl/rgm) to generate random geological models. 

The code is written in Fortran + MPI. Currently, it can only be compiled with Intel's compiler suite [Intel oneAPI HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html).

# Use
To install `LATTE`, 

```
cd src
ruby install.rb
```

The compiled `LATTE` executables will be at `bin`.

To remake, 

```
cd src
ruby install.rb clean
```

We include several simple examples to use `LATTE` in [example](example). To run the tests, 

```
cd test
```

and the scripts to reproduce the examples in the mansucript are contained in subfolders. 

# Examples
Reproducible examples associated with [the paper](https://academic.oup.com/gji/article/241/2/1275/8046728) is in the `example` directory. 

# License
&copy; 2024-2026. Triad National Security, LLC. All rights reserved. 

This program is Open-Source under the BSD-3 License.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
- Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Author
Kai Gao, <kaigao@lanl.gov>

We welcome feedback, bug reports, and suggestions for improving `LATTE`. 

If you use this package in your research and find it useful, please cite it as

* Kai Gao, Ting Chen, 2025, LATTE: open-source, high-performance traveltime computation, tomography and source location in acoustic and elastic media, Geophysical Journal International, doi: [10.1093/gji/ggaf079](https://academic.oup.com/gji/article/241/2/1275/8046728). 
* Kai Gao, Ting Chen, 2024, LATTE: Los Alamos TravelTime package based on Eikonal equation, GitHub Repository, url: [github.com/lanl/latte_traveltime](https://github.com/lanl/latte_traveltime)

