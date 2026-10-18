#!/bin/bash

# Compare x_eikonal3, x_fatt3, x_trtt3 and x_tloc3 with their CUDA versions:
# runtime and byte-by-byte identity of the results. Needs a build with
# use_cuda = on and a GPU. bindir and np can be set in the environment.

set -e

bindir=${bindir:-../../bin}
np=${np:-4}

python3 compare.py --bindir $bindir --np $np
