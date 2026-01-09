
# paths
bindir = $(HOME)/bin
objdir = $(bindir)/obj/select
moddir = $(bindir)/obj/select

# dependencies
object = main_select.o

obj = $(addprefix $(objdir)/,$(object))

# targets
exec = $(bindir)/x_select

all: makedir $(exec)

# options
flitdir = $(HOME)/src/flit

include $(flitdir)/src/Makefile.in
inc = $(base_inc) -I$(flitdir)/lib
lflags = $(flitdir)/lib/libflit.a $(base_lflags)
fflags = $(base_fflags)
ccflags = $(base_ccflags)

# compile
$(objdir)/%.o : ./%.f90
	$(fc) -o $@ -c $(fflags) $(inc) $<

# link
$(exec) : $(obj)
	$(fc) -o $@ $^ $(lflags) $(inc)

# make directory
makedir:
	-@mkdir -p $(bindir)
	-@mkdir -p $(objdir)
	-@mkdir -p $(moddir)

# clean
clean:
	-@rm -rf $(objdir)/*.o $(moddir)/*.mod
