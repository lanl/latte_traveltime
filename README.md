# Description
**LATTE: Los Alamos TravelTime package based on Eikonal equation**

We develop an open-source package for traveltime computation, tomography and source location in 2D/3D isotropic acoustic/elastic media. 

The work was supported by Los Alamos National Laboratory (LANL) Laboratory Directory Research and Development (LDRD) project 20240322ER. LANL is operated by Triad National Security, LLC, for the National Nuclear Security Administration (NNSA) of the U.S. Department of Energy (DOE) under Contract No. 89233218CNA000001. The research used high-performance computing resources provided by LANL's Institutional Computing program. 

The work is under LANL open source approval reference O4770.

# Reference
Please refer to our submitted manuscript [LA-UR-24-27806](doc/manuscript_submitted_gji.pdf) for details. 

A mannual will be released in the near future pending separate approval.

# Requirement
`LATTE` depends on [FLIT](https://github.com/lanl/flit). Some examples in [example](example) use [RGM](https://github.com/lanl/rgm) to generate random geological models. 

The code is written in Fortran. Currently, it only support Intel's Fortran compiler, which is freely available through [Intel oneAPI Base Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit.html#gs.bed72v) and [Intel HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html#gs.bed5op). 

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

# License
&copy; 2024. Triad National Security, LLC. All rights reserved. 

This program is Open-Source under the BSD-3 License.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
- Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Author
Kai Gao, <kaigao@lanl.gov>

We welcome feedbacks, bug reports, and improvement ideas on `LATTE`. 

If you use this package in your research and find it useful, please cite it as

* Kai Gao, Ting Chen, 2024, LATTE: Los Alamos TravelTime package based on Eikonal equation, url: [github.com/lanl/latte](https://github.com/lanl/latte)
* Kai Gao, Ting Chen, 2024, LATTE: Open-source, high-performance acoustic and elastic traveltime computation, tomography, and source location, _under review_. 
