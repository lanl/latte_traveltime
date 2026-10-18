#
# © 2024-2026. Triad National Security, LLC. All rights reserved.
#
# This program was produced under U.S. Government contract 89233218CNA000001
# for Los Alamos National Laboratory (LANL), which is operated by
# Triad National Security, LLC for the U.S. Department of Energy/National Nuclear
# Security Administration. All rights in the program are reserved by
# Triad National Security, LLC, and the U.S. Department of Energy/National
# Nuclear Security Administration. The Government is granted for itself and
# others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
# license in this material to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display publicly,
# and to permit others to do so.
#
# Author:
#    Kai Gao, kaigao@lanl.gov
#

# paths
bindir = $(PWD)/../../bin
objdir = $(HOME)/bin/obj/latte3
moddir = $(HOME)/bin/obj/latte3

# Objects and modules of the CUDA executables, compiled with -Dcuda
objdir_cuda := $(objdir)_cuda
moddir_cuda := $(moddir)_cuda

# Dependencies
common = \
    module_parameters.o \
	module_vars.o \
	module_utility.o \
	module_traveltime_iso.o \
	module_traveltime_iso_reflection.o

# The CUDA executables add the CUDA solvers, which module_traveltime_iso uses
common_cuda = \
	module_parameters.o \
	module_vars.o \
	module_utility.o \
	traveltime_iso_cuda.o \
	module_traveltime_iso_cuda.o \
	module_traveltime_iso.o \
	module_traveltime_iso_reflection.o

object1 = \
    $(common) \
	main_eikonal.o

common_inversion = \
    module_gradient.o \
    submodule_gradient_fatt.o \
    submodule_gradient_trtt.o \
    submodule_gradient_tloc.o \
    submodule_reflector_imaging.o \
	module_regularization.o \
	module_inversion_regularization.o \
	module_inversion_gradient.o \
	module_inversion_search_direction.o \
	module_inversion_step_size.o

object2 = \
	$(common) \
    $(common_inversion) \
	main_fatt.o

object3 = \
	$(common) \
    $(common_inversion) \
	main_trtt.o

object4 = \
	$(common) \
    $(common_inversion) \
	main_tloc.o

obj1 = $(addprefix $(objdir)/, $(object1))
obj2 = $(addprefix $(objdir)/, $(object2))
obj3 = $(addprefix $(objdir)/, $(object3))
obj4 = $(addprefix $(objdir)/, $(object4))

obj1_cuda = $(addprefix $(objdir_cuda)/, $(common_cuda) main_eikonal.o)
obj2_cuda = $(addprefix $(objdir_cuda)/, $(common_cuda) $(common_inversion) main_fatt.o)
obj3_cuda = $(addprefix $(objdir_cuda)/, $(common_cuda) $(common_inversion) main_trtt.o)
obj4_cuda = $(addprefix $(objdir_cuda)/, $(common_cuda) $(common_inversion) main_tloc.o)

# targets
exec1 = $(bindir)/x_eikonal3
exec2 = $(bindir)/x_fatt3
exec3 = $(bindir)/x_trtt3
exec4 = $(bindir)/x_tloc3

exec1_cuda = $(bindir)/x_eikonal3_cuda
exec2_cuda = $(bindir)/x_fatt3_cuda
exec3_cuda = $(bindir)/x_trtt3_cuda
exec4_cuda = $(bindir)/x_tloc3_cuda

all: makedir $(exec1) $(exec2) $(exec3) $(exec4)

# options
flitdir = $(HOME)/src/flit

include $(flitdir)/src/Makefile.in
inc = $(base_inc) -I$(flitdir)/lib
lflags = $(flitdir)/lib/libflit.a $(base_lflags)
fflags = $(base_fflags) -Ddim3 \
	#	-O0 -g -check bounds

# CUDA executables, off by default. With use_cuda = on here or make
# use_cuda=on, the build adds x_eikonal3_cuda, x_fatt3_cuda, x_trtt3_cuda and
# x_tloc3_cuda, which solve the eikonal and adjoint-state equations on GPUs,
# to the CPU executables, which stay the same. cuda_arch is the GPU
# architecture, e.g., sm_80 or sm_90 (nvcc --list-gpu-arch), and native
# targets the GPU of the build machine; run make clean after a change. Set
# cuda_ccbin to a host compiler that nvcc supports, e.g., g++-12, if the
# default one is not. cuda_home is the CUDA Toolkit directory whose lib64
# holds the CUDA runtime; the default is the parent of the directory of nvcc.
use_cuda ?= off
cuda_arch ?= native
cuda_ccbin ?=
nvcc = nvcc
cuda_home ?= $(abspath $(dir $(shell which $(nvcc)))..)

ifeq ($(use_cuda), on)
all: makedir_cuda $(exec1_cuda) $(exec2_cuda) $(exec3_cuda) $(exec4_cuda)
endif

lflags_cuda = $(lflags) -L$(cuda_home)/lib64 -lcudart_static -lrt -ldl -lpthread
nvccflags = -O3 --fmad=false -arch=$(cuda_arch) -Xcompiler -fPIC \
	$(if $(cuda_ccbin), -ccbin $(cuda_ccbin))

# The CUDA objects define cuda and keep their modules apart from the CPU ones,
# also when moddir or fflags come from the command line
$(objdir_cuda)/%.o : override moddir := $(moddir_cuda)
$(objdir_cuda)/%.o : override fflags += -Dcuda

# compile
$(objdir)/%.o : $(PWD)/%.f90
	$(fc) -o $@ -c $(fflags) $(inc) $<

$(objdir)/%.o : $(PWD)/../%.f90
	$(fc) -o $@ -c $(fflags) $(inc) $<

$(objdir_cuda)/%.o : $(PWD)/%.f90
	$(fc) -o $@ -c $(fflags) $(inc) $<

$(objdir_cuda)/%.o : $(PWD)/../%.f90
	$(fc) -o $@ -c $(fflags) $(inc) $<

# --fmad=false keeps the GPU arithmetic the same as the CPU arithmetic
$(objdir_cuda)/%.o : $(PWD)/%.cu
	$(nvcc) -o $@ -c $(nvccflags) $<

# link
$(exec1) : $(obj1)
	$(fc) -o $@ $^ $(lflags) $(inc)

$(exec2) : $(obj2)
	$(fc) -o $@ $^ $(lflags) $(inc)

$(exec3) : $(obj3)
	$(fc) -o $@ $^ $(lflags) $(inc)

$(exec4) : $(obj4)
	$(fc) -o $@ $^ $(lflags) $(inc)

$(exec1_cuda) : $(obj1_cuda)
	$(fc) -o $@ $^ $(lflags_cuda) $(inc)

$(exec2_cuda) : $(obj2_cuda)
	$(fc) -o $@ $^ $(lflags_cuda) $(inc)

$(exec3_cuda) : $(obj3_cuda)
	$(fc) -o $@ $^ $(lflags_cuda) $(inc)

$(exec4_cuda) : $(obj4_cuda)
	$(fc) -o $@ $^ $(lflags_cuda) $(inc)

# make directory
makedir:
	-@mkdir -p $(bindir)
	-@mkdir -p $(objdir)
	-@mkdir -p $(moddir)

makedir_cuda:
	-@mkdir -p $(objdir_cuda)
	-@mkdir -p $(moddir_cuda)

# clean
clean:
	-@rm -rf $(objdir)/*.o $(moddir)/*.mod $(moddir)/*.smod \
		$(objdir_cuda)/*.o $(moddir_cuda)/*.mod $(moddir_cuda)/*.smod \
		$(exec1) $(exec2) $(exec3) $(exec4) \
		$(exec1_cuda) $(exec2_cuda) $(exec3_cuda) $(exec4_cuda)

