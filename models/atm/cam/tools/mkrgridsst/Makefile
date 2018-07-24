EXENAME = mkrgridsst
FC = f90
RM = rm
VPATH = /fs/cgd/home0/rosinski/mkrgridsst
RM := rm

# Determine platform 
UNAMES := $(shell uname -s)
UNAMEM := $(findstring CRAY,$(shell uname -m))

# Platform-independent settings

EXEDIR = .
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf

# Check for NetCDF library and include directories 

ifeq ($(LIB_NETCDF),$(null))
LIB_NETCDF := /usr/local/lib
endif

ifeq ($(INC_NETCDF),$(null))
INC_NETCDF := /usr/local/include
endif


# Architecture-specific flags and rules
#
#------------------------------------------------------------------------
# Cray 
#------------------------------------------------------------------------

ifeq ($(UNAMEM),CRAY)
endif

#------------------------------------------------------------------------
# SGI
#------------------------------------------------------------------------

ifeq ($(UNAMES),IRIX64)
CFLAGS = -64 -DSGI -g
FFLAGS = -64 -DSGI -g -c -cpp -trapuv -C -I$(INC_NETCDF)
LDFLAGS += -64 
endif

#------------------------------------------------------------------------
# SUN
#------------------------------------------------------------------------

ifeq ($(UNAMES),SunOS)
CFLAGS = -DSUN -g
FFLAGS = -DSUN -xpp=cpp -g -c -stackvar -f -I$(INC_NETCDF)
endif

#------------------------------------------------------------------------
# Linux
#------------------------------------------------------------------------

ifeq ($(UNAMES),Linux)
FC = pgf90
FFLAGS = -fast -g -c -I$(INC_NETCDF)
endif

#------------------------------------------------------------------------
# Default rules and macros
#------------------------------------------------------------------------

OBJS := chkdims.o fmain.o lininterp.o wrap_nf.o

.SUFFIXES:
.SUFFIXES: .F90 .o

.F90.o:
	$(FC) $(FFLAGS) $<


$(EXEDIR)/$(EXENAME): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

junk: fourier.o junk.o
	$(FC) -o $@ fourier.o junk.o $(LDFLAGS)

clean:
	$(RM) -f $(OBJS) $(EXEDIR)/$(EXENAME)

chkdims.o: chkdims.F90
fmain.o: fmain.F90
lininterp.o: lininterp.F90
wrap_nf.o: wrap_nf.F90
