# Makefile to build mkrgrid on various platforms
# Note: If netcdf library is not built in the standard location, you must set the environment
# variables INC_NETCDF and LIB_NETCDF
#
# To get accuracy to 64-bits from FFT, need to specify real*8 shr_kind_mod in
# constants.  Thus FFLAGS_FIXED contains an argument to assure this.  FFLAGS
# does not need this option because the code is written with explicit shr_kind_mod
# defined for all floating point variables.

EXEDIR = .
EXENAME = mkrgrid
RM = rm

.SUFFIXES:
.SUFFIXES: .f90 .F .o

# Check for the NetCDF library and include directories 
ifeq ($(LIB_NETCDF),$(null))
LIB_NETCDF := /usr/local/lib
endif

ifeq ($(INC_NETCDF),$(null))
INC_NETCDF := /usr/local/include
endif

# Determine platform 
UNAMES := $(shell uname -s)
UNAMEM := $(findstring CRAY,$(shell uname -m))

# Architecture-specific flags and rules
#
#------------------------------------------------------------------------
# Cray 
#------------------------------------------------------------------------

ifeq ($(UNAMEM),CRAY)
FC = f90
FFLAGS = -c -I$(INC_NETCDF)
FFLAGS_FIXED = -c -f fixed
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# SGI
#------------------------------------------------------------------------

ifeq ($(UNAMES),IRIX64)
FC = f90
FFLAGS = -64 -O2 -c -trapuv -I$(INC_NETCDF)
FFLAGS_FIXED = -64 -r8 -O2 -c -trapuv -fixedform
LDFLAGS = -64 -L/usr/local/lib64/r4i4 -lnetcdf
endif

#------------------------------------------------------------------------
# SUN
#------------------------------------------------------------------------

ifeq ($(UNAMES),SunOS)
FC = f90
FFLAGS = -c -stackvar -f -I$(INC_NETCDF)
FFLAGS_FIXED = -c -stackvar -f -fixed
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# AIX
#------------------------------------------------------------------------

ifeq ($(UNAMES),AIX)
FC = xlf90
FFLAGS = -c -qsuffix=f=f90 -qfree=f90 -I$(INC_NETCDF)
FFLAGS_FIXED = -c -qfixed -qdpc
LDFLAGS = -L/usr/local/lib32/r4i4 -lnetcdf
endif

#------------------------------------------------------------------------
# OSF1
#------------------------------------------------------------------------

ifeq ($(UNAMES),OSF1)
FC = f90
FFLAGS = -c -I$(INC_NETCDF)
FFLAGS_FIXED = -c -fixed -r8
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# Linux
#------------------------------------------------------------------------

ifeq ($(UNAMES),Linux)
FC = pgf90
FFLAGS =       -c -I$(INC_NETCDF) -fast
FFLAGS_FIXED = -c -Mnofreeform -fast
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# Default rules and macros
#------------------------------------------------------------------------

OBJS := shr_kind_mod.o constants.o fmain.o mkrgrid.o rgconvert.o fourier.o \
	lininterp.o cubinterp.o fortfft.o wrap_nf.o control.o gridspecs.o \
        interp_driver.o timing_stubs.o

.SUFFIXES:
.SUFFIXES: .f90 .F .o

.f90.o:
	$(FC) $(FFLAGS) $<

.F.o:
	$(FC) $(FFLAGS_FIXED) $<

$(EXEDIR)/$(EXENAME): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

clean:
	$(RM) -f $(OBJS) *.mod $(EXEDIR)/$(EXENAME)

shr_kind_mod.o: shr_kind_mod.f90
constants.o: constants.f90 shr_kind_mod.o
control.o: control.f90 constants.o
gridspecs.o: gridspecs.f90 constants.o
fourier.o: fourier.f90 fortfft.o shr_kind_mod.o constants.o
lininterp.o: lininterp.f90 shr_kind_mod.o constants.o
fmain.o: fmain.f90 constants.o control.o
mkrgrid.o: mkrgrid.f90 constants.o control.o gridspecs.o
rgconvert.o: rgconvert.f90 shr_kind_mod.o control.o
cubinterp.o: cubinterp.f90 shr_kind_mod.o constants.o
wrap_nf.o: wrap_nf.f90 shr_kind_mod.o
interp_driver.o: interp_driver.f90 shr_kind_mod.o control.o gridspecs.o
fortfft.o: fortfft.F shr_kind_mod.o
timing_stubs.o: timing_stubs.f90
