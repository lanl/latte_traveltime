#!/bin/bash

# Test the level ranges of the 3D fast sweeps on a grid with ny < nz. bindir
# and np can be set in the environment.

set -e

bindir=${bindir:-../../bin}
np=${np:-3}

python3 test.py --bindir $bindir --np $np
