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
# license in this material to reproduce, prepare. derivative works, 
# distribute copies to the public, perform publicly and display publicly, 
# and to permit others to do so.
#
# Author:
#    Kai Gao, kaigao@lanl.gov
#

# paths
bindir = $(PWD)
objdir = ./obj
moddir = ./obj

# dependencies
object2 = generate_data_2d.o
object3 = generate_data_3d.o

obj2 = $(addprefix $(objdir)/, $(object2))
obj3 = $(addprefix $(objdir)/, $(object3))

# targets
exec2 = $(PWD)/exec2
exec3 = $(PWD)/exec3

all: makedir $(exec2) $(exec3)

# options
flitdir = $(HOME)/src/libflit
rgmdir = $(HOME)/src/rgm

include $(flitdir)/src/Makefile.in
inc = $(base_inc) \
	-I$(flitdir)/lib \
	-I$(rgmdir)/lib
lflags = $(rgmdir)/lib/librgm.a \
	$(flitdir)/lib/libflit.a \
	$(base_lflags)
fflags = $(base_fflags) #-O0 -g -check bounds

# compile
$(objdir)/%.o : ./%.f90
	$(fc) -o $@ -c $(fflags) $(inc) $<

# link
$(exec2) : $(obj2)
	$(fc) -o $@ $^ $(lflags) $(inc)

$(exec3) : $(obj3)
	$(fc) -o $@ $^ $(lflags) $(inc)

# make directory
makedir:
	-@mkdir -p $(bindir)
	-@mkdir -p $(objdir)
	-@mkdir -p $(moddir)

# clean
clean:
	-@rm -rf $(objdir)/*.o $(moddir)/*.mod
