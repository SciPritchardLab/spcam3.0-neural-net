# Makefile to build interpic on various platforms
# Note: If netcdf library is not built in the standard location, you must set the environment
# variables INC_NETCDF and LIB_NETCDF

# Set up special characters
null  :=

EXEDIR = .
EXENAME = interpaerosols
RM = rm

.SUFFIXES:
.SUFFIXES: .f90 .o

# Check for the NetCDF library and include directories 
ifeq ($(LIB_NETCDF),$(null))
LIB_NETCDF := /usr/local/lib
endif

ifeq ($(INC_NETCDF),$(null))
INC_NETCDF := /usr/local/include
endif

# Determine platform 
UNAMES := $(shell uname -s)

# Architecture-specific flags and rules
#------------------------------------------------------------------------
# SGI
#------------------------------------------------------------------------

ifeq ($(UNAMES),IRIX64)
FC      = f90
FFLAGS  = -64 -c -r8 -i4 -I$(INC_NETCDF) -O2
LDFLAGS = -64 -L$(LIB_NETCDF) -lnetcdf -lscs
endif

#------------------------------------------------------------------------
# AIX
#------------------------------------------------------------------------

ifeq ($(UNAMES),AIX)
FC      = xlf90
FFLAGS  = -c -I$(INC_NETCDF) -q64 -qsuffix=f=f90 -O2 -qmaxmem=-1
LDFLAGS = -L$(LIB_NETCDF) -q64 -lnetcdf -lblas
#LDFLAGS = -L$(LIB_NETCDF) -q64 -lnetcdf -llapack -lblas
endif

#------------------------------------------------------------------------
# OSF1
#------------------------------------------------------------------------

ifeq ($(UNAMES),OSF1)
FC      = f90
FFLAGS  = -c -I$(INC_NETCDF)
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf -lcxml
endif

#------------------------------------------------------------------------
# Linux
#------------------------------------------------------------------------

ifeq ($(UNAMES),Linux)
  ifeq ($(USER_FC),$(null))
    FC := pgf90
  else
    FC := $(USER_FC)
  endif

  ifeq ($(FC),pgf90)
    ifeq ($(DEBUG),TRUE)
      FFLAGS = -c -I$(INC_NETCDF) -g -Ktrap=fp -Mbounds
    else
      FFLAGS = -c -I$(INC_NETCDF) -fast
    endif
    LDFLAGS = -L$(LIB_NETCDF) -lnetcdf -llapack -lblas
  endif

  ifeq ($(FC),lf95)
    ifeq ($(DEBUG),TRUE)
      FFLAGS =  -c -I$(INC_NETCDF) -g --chk e,s,u
    else
      FFLAGS =  -c -I$(INC_NETCDF) -O
    endif
    LDFLAGS = -L$(LIB_NETCDF) -lnetcdf -llapackmt -lblasmt
  endif
endif

#------------------------------------------------------------------------
# Default rules and macros
#------------------------------------------------------------------------

OBJS := addglobal.o bilin.o driver.o fmain.o globals.o prec.o preserve_mean.o wrap_nf.o

.f90.o:
	$(FC) $(FFLAGS) $<

$(EXEDIR)/$(EXENAME): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

clean:
	$(RM) -f $(OBJS) *.mod $(EXEDIR)/$(EXENAME)

addglobal.o: 
bilin.o: prec.o
driver.o: prec.o globals.o preserve_mean.o
fmain.o: globals.o
globals.o:
prec.o:
preserve_mean.o: prec.o globals.o
wrap_nf.o: prec.o

