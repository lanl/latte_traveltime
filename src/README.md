
We currently only provide traditional Makefile-based approach to compile `LATTE`:

* To compile, set the path to `FLIT` with `flitdir = ...` in [Makefile](2d/Makefile) and [Makefile](3d/Makefile); the default path to `FLIT` is `$HOME/src/flit`. 
* To use legacy eikonal solvers, enable `-Dlegacy_solver` in [Makefile](2d/Makefile) and [Makefile](3d/Makefile); the default is not using legacy solvers. 
