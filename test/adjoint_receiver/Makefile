
# Test drivers that call forward_iso and adjoint_iso of LATTE directly. They
# link the objects and modules of a LATTE build: latte2 and latte3 are the
# object directories of src/2d and src/3d (objdir in their Makefiles), which
# hold the CPU objects also when the build adds the CUDA executables.

flitdir = $(HOME)/src/flit
latte2 = $(HOME)/bin/obj/latte2
latte3 = $(HOME)/bin/obj/latte3

moddir = .

include $(flitdir)/src/Makefile.in
inc = $(base_inc) -I$(flitdir)/lib
lflags = $(flitdir)/lib/libflit.a $(base_lflags)

objects = module_parameters.o module_vars.o module_utility.o module_traveltime_iso.o

all: kernel2d kernel3d

kernel2d: kernel2d.f90
	$(fc) -o $@ $(base_fflags) -Ddim2 -I$(latte2) $(inc) $< $(addprefix $(latte2)/, $(objects)) $(lflags)

kernel3d: kernel3d.f90
	$(fc) -o $@ $(base_fflags) -Ddim3 -I$(latte3) $(inc) $< $(addprefix $(latte3)/, $(objects)) $(lflags)

clean:
	-@rm -f kernel2d kernel3d *.mod

