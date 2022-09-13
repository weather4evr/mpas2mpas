SHELL           =       /bin/sh
DEVTOP          =       `pwd`
LIBINCLUDE      =       .
.SUFFIXES: .F .i .o .f90 .c

#### Architecture specific settings ####

F90             = ifort # mpif90
LD              = ifort #mpif90
F77             = ifort
CC              = gcc
CPP             = cpp -P

FCFFLAGS        = -g -assume byterecl -heap-arrays #-convert big_endian
CCFLAGS         = -O3
CPPFLAGS        =
LDFLAGS         = $(FCFFLAGS)
NETCDFPATH      = $(NETCDF)
# PNETCDF needed to output CDF-5 variables
PNETCDFPATH     = $(PNETCDF)
OPTIMIZATION    = -O3 -xHOST
DEBUG           = -traceback #-C

MODS = kinds.o                             \
       namelist_mod.o                      \
       mpas_netcdf_interface.o             \
       interpolation_stuff.o               

SUBS =

# Define compile time variables

NETCDFcomp = -L$(NETCDFPATH)/lib -lnetcdf -lnetcdff -I$(NETCDFPATH)/include -L$(PNETCDFPATH)/lib -lpnetcdf -I$(PNETCDFPATH)/include

# Define executables

all: mpas2mpas

mpas2mpas:           mpas2mpas.f90  $(MODS) $(SUBS)
	$(LD) $(LDFLAGS)           mpas2mpas.f90  $(MODS) $(SUBS) $(NETCDFcomp) $(OPTIMIZATION) $(DEBUG) -lm -o $@.exe

clean:
	-rm *.o *~ \#* *.mod *.bin fort.* *.exe >& /dev/null

# Define make rules

.f90.o:
	$(F90) $(CPPFLAGS) $(FCFFLAGS) $(NETCDFcomp) $(OPTIMIZATION) $(DEBUG) -c $*.f90

.f.o:
	$(F77) $(CPPFLAGS) $(FCFFLAGS) $(NETCDFcomp) $(OPTIMIZATION) $(DEBUG) -c $*.f

.F.o:
	$(CPP) $(CPPFLAGS) $< > $*.f90
	$(F90) $(CPPFLAGS) $(FCFFLAGS) $(NETCDFcomp) $(OPTIMIZATION) $(DEBUG) -c $*.f90

