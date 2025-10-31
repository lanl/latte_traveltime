#
# © 2024. Triad National Security, LLC. All rights reserved.
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
objdir = $(HOME)/bin/obj/latte2
moddir = $(HOME)/bin/obj/latte2

# Dependencies
common = \
    module_parameters.o \
	module_vars.o \
	module_utility.o \
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

# targets
exec1 = $(bindir)/x_eikonal2
exec2 = $(bindir)/x_fatt2
exec3 = $(bindir)/x_trtt2
exec4 = $(bindir)/x_tloc2

all: makedir $(exec1) $(exec2) $(exec3) $(exec4)

# options
flitdir = $(HOME)/src/flit

include $(flitdir)/src/Makefile.in
inc = $(base_inc) -I$(flitdir)/lib
lflags = $(flitdir)/lib/libflit.a $(base_lflags)
fflags = $(base_fflags) -Ddim2 \
	# -O0 -g -check bounds
	# -Dlegacy_solver

# compile
$(objdir)/%.o : $(PWD)/%.f90
	$(fc) -o $@ -c $(fflags) $(inc) $<

$(objdir)/%.o : $(PWD)/../%.f90
	$(fc) -o $@ -c $(fflags) $(inc) $<

# link
$(exec1) : $(obj1)
	$(fc) -o $@ $^ $(lflags) $(inc)

$(exec2) : $(obj2)
	$(fc) -o $@ $^ $(lflags) $(inc)

$(exec3) : $(obj3)
	$(fc) -o $@ $^ $(lflags) $(inc)

$(exec4) : $(obj4)
	$(fc) -o $@ $^ $(lflags) $(inc)

# make directory
makedir:
	-@mkdir -p $(bindir)
	-@mkdir -p $(objdir)
	-@mkdir -p $(moddir)

# clean
clean:
	-@rm -rf $(objdir)/*.o $(moddir)/*.mod $(moddir)/*.smod \
        $(exec1) $(exec2) $(exec3) $(exec4)

