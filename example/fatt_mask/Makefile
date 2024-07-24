
# paths
bindir = $(PWD)
objdir = ./obj
moddir = ./obj

# dependencies
object = create_model_and_geometry.o

obj = $(addprefix $(objdir)/, $(object))

# targets
exec = $(PWD)/exec

all: makedir $(exec)

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
