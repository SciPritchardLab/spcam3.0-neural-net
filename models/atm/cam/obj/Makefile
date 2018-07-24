# Make macros for CAM.

ROOTDIR      := /Users/mark/mmf/cam3_sp
EXENAME      := cam
MODEL_EXEDIR := /Users/mark/mmf/cam3_sp/models/atm/cam/run
INC_NETCDF   := /usr/local/absoft/include
LIB_NETCDF   := /usr/local/absoft/lib
MOD_NETCDF   := 
INC_MPI      := 
LIB_MPI      := 
ESMF_ROOT    := /Users/mark/mmf/cam3_sp/models/utils/esmf
ESMF_BLD     := /Users/mark/mmf/cam3_sp/models/atm/cam/obj/esmf
DEBUG        := FALSE
SMP          := FALSE
USER_FC      := 
USER_CC      := 
USER_CFLAGS  := 
USER_FFLAGS  := 
USER_LDFLAGS := 

#-----------------------------------------------------------------------
# Makefile for building CAM on various platforms.
#
# This makefile assumes the existence the file:
#
# Filepath      The directories, listed one per line, that contain the source
#               code required to build CAM.  This list is used to set the
#               VPATH variable which is used by GNU make to search for
#               dependencies (after it looks in the directory from which
#               it was invoked).  This list of directories, prepended with ./,
#               is also used to construct the list of search directories 
#               used by the preprocessor (as specified by -I command-line options).
#
# The following macros may be set in the user's environment:
#

# ROOTDIR       Root directory for the CAM source distribution.  If ROOTDIR is not
#               set the makefile searches for it in a file called "Rootdir".
# EXENAME       Name to call the executable.  Default: atm
# MODEL_EXEDIR  Directory to build the executable in.  Default: ./
# INC_NETCDF    Directory containing the NetCDF include files.  Default: /usr/local/include
# LIB_NETCDF    Directory containing the NetCDF library.  Default: /usr/local/lib
# MOD_NETCDF    Directory containing the NetCDF module files.  Default: /usr/local/lib
# INC_MPI       Directory containing the MPI include files.  Default: /usr/local/include
# LIB_MPI       Directory containing the MPI library.  Default: /usr/local/lib
# ESMF_ROOT     Root directory of ESMF source tree.  Default: $(ROOTDIR)/models/utils/esmf
# ESMF_BLD      Upper level build directory for ESMF
#               The ESMF .o files are built in $(ESMF_BLD)/obj
#               The ESMF .m files are built in arch dependent subdirs of $(ESMF_BLD)/mod
#               The ESMF library is built in architecture and optimization
#               dependent subdirectories of $(ESMF_BLD)/lib
#               Default: ./esmf
# DEBUG         Set to TRUE to turn on compiler debugging options.  Default: FALSE 
# SPMD          Whether to build in SPMD mode or not. 	[values TRUE FALSE]
#               Default: read from ./misc.h
#               N.B. To ensure consistency between the Makefile and misc.h the SPMD
#                    macro should not be set.  This forces it to be read from misc.h.
# SMP           Set to TRUE to enable building in SMP mode (uses OpenMP).  
#               Currently implemented for IBM, SGI, linux-pgf90. (default is TRUE on IBM and 
#               linux-pgf90, and depends on SPMD setting on SGI).
# USER_FC       Allow user to override the default Fortran compiler specified in Makefile.
# USER_CC       Allow user to override the default C compiler specified in Makefile (linux only).
# USER_CFLAGS   Additional C compiler flags that the user wishes to set.
# USER_FFLAGS   Additional Fortran compiler flags that the user wishes to set.
# USER_LDLAGS   Additional load flags that the user wishes to set.
# F_OPTIMIZATION_OVERRIDE
#               Use this to replace default Fortran compiler optimization.  
# NO_SWITCH     On Compaq if the hardware switch is not available
#               set this env variable to "TRUE".
#
# Note: The ESMF library is included in the CAM distribution in 
#       $ROOTDIR/models/utils/esmf and is built using this makefile.  
#------------------------------------------------------------------------

# Set up special characters
null  :=
space := $(null) $(null)
comma := $(null),$(null)

# Determine distribution root directory
ifeq ($(ROOTDIR),$(null))
ROOTDIR := $(shell cat Rootdir)
endif

# Check for override of default Fortran compiler optimizations
ifneq ($(F_OPTIMIZATION_OVERRIDE),$(null))
FORTRAN_OPTIMIZATION := $(F_OPTIMIZATION_OVERRIDE)
endif

# Check for the NetCDF library and include directories 
ifeq ($(LIB_NETCDF),$(null))
LIB_NETCDF := /usr/local/lib
endif

ifeq ($(INC_NETCDF),$(null))
INC_NETCDF := /usr/local/include
endif

ifeq ($(MOD_NETCDF),$(null))
MOD_NETCDF := $(LIB_NETCDF)
endif

# Check for the MPI library and include directories 
ifeq ($(LIB_MPI),$(null))
LIB_MPI := /usr/local/lib
endif

ifeq ($(INC_MPI),$(null))
INC_MPI := /usr/local/include
endif

# Build the ESMF library
ifeq ($(ESMF_ROOT),$(null))
ESMF_ROOT := ${ROOTDIR}/models/utils/esmf
endif

ifeq ($(ESMF_BLD),$(null))
ESMF_BLD := $(shell (mkdir ./esmf 2>/dev/null;cd ./esmf;pwd))
endif

ESMF_BOPT := O
ifeq ($(DEBUG),TRUE)
  ESMF_BOPT := g
endif
ESMF_MOD := $(ESMF_BLD)/mod/mod$(ESMF_BOPT)
ESMF_LIB := $(ESMF_BLD)/lib/lib$(ESMF_BOPT)

# Check for directory in which to put executable
ifeq ($(MODEL_EXEDIR),$(null))
MODEL_EXEDIR := .
endif

# Check for name of executable
ifeq ($(EXENAME),$(null))
EXENAME := atm
endif

# Check if SPMD is defined in "misc.h"
# Ensure that it is defined and not just "undef SPMD" set in file
ifeq ($(SPMD),$(null))
	SPMDSET := $(shell /bin/grep SPMD misc.h)
	ifneq (,$(findstring define,$(SPMDSET)))
		SPMD := TRUE
	else
		SPMD := FALSE
	endif
endif

# Load dependency search path.
dirs := . $(shell cat Filepath)

# Set cpp search path, include netcdf
cpp_dirs := $(dirs) $(INC_NETCDF) $(INC_MPI)
cpp_path := $(foreach dir,$(cpp_dirs),-I$(dir)) # format for command line

# Expand any tildes in directory names. Change spaces to colons.
VPATH    := $(foreach dir,$(cpp_dirs),$(wildcard $(dir))) 
VPATH    := $(subst $(space),:,$(VPATH))               

#------------------------------------------------------------------------
# Primary target: build the model
#------------------------------------------------------------------------
all: $(MODEL_EXEDIR)/$(EXENAME)

# Get list of files and build dependency file for all .o files
#   using perl scripts mkSrcfiles and mkDepends

SOURCES   := $(shell cat Srcfiles)

# Newer makes set the CURDIR variable.
CURDIR := $(shell pwd)

$(CURDIR)/Depends: $(CURDIR)/Srcfiles $(CURDIR)/Filepath
	$(ROOTDIR)/models/atm/cam/bld/mkDepends Filepath Srcfiles > $@

$(CURDIR)/Srcfiles: $(CURDIR)/Filepath
	$(ROOTDIR)/models/atm/cam/bld/mkSrcfiles > $@

OBJS      := $(addsuffix .o, $(basename $(SOURCES)))

CPPDEF    := -DCAM -DNO_SHR_VMATH -DHIDE_SHR_MSG

$(MODEL_EXEDIR)/$(EXENAME): $(OBJS)
	$(FC) -o $@ $(OBJS) -L$(LIB_NETCDF) -lnetcdf -L$(ESMF_LIB)/$(ESMF_ARCH) -lesmf $(LDFLAGS)

debug: $(OBJS)
	echo "FFLAGS: $(FFLAGS)"
	echo "LDFLAGS: $(LDFLAGS)"
	echo "OBJS: $(OBJS)"

test_fc: test_fc.o
	$(FC) -o $@ test_fc.o $(LDFLAGS)
test_nc: test_nc.o
	$(FC) -o $@ test_nc.o -L$(LIB_NETCDF) -lnetcdf $(LDFLAGS)
test_mpi: test_mpi.o
	$(FC) -o $@ test_mpi.o $(LDFLAGS)

# Architecture-specific flags and rules
#
# Determine platform 
UNAMES := $(shell uname -s)

#------------------------------------------------------------------------
# AIX
#------------------------------------------------------------------------

ifeq ($(UNAMES),AIX)

ifeq ($(SMP),$(null))
  SMP := TRUE
endif

#ADDRESS    := Q32
 ADDRESS    := Q64
ifeq ($(ADDRESS),Q32)
  QCMP      :=
  QLDR      := -bmaxdata:0x80000000 -bmaxstack:0x10000000
  ESMF_ARCH := rs6000_sp
endif
ifeq ($(ADDRESS),Q64)
  QCMP      := -q64
  QLDR      := -q64
  ESMF_ARCH := rs6000_64
endif

CPP        := /lib/cpp
CC         := xlc_r
CFLAGS     := $(cpp_path) -O2 -DAIX $(CPPDEF) $(QCMP)
FPPFLAGS   := -WF,-DAIX,$(subst $(space),$(comma),$(CPPDEF))
#
# Activate the line below for MPI derived types in FV transposes (AAM)
# FPPFLAGS   += -WF,-DUSE_MPI_TYPES
#
mod_path   := -I$(ESMF_MOD)/$(ESMF_ARCH) -I$(MOD_NETCDF)
FFLAGS     := $(cpp_path) $(mod_path) $(FPPFLAGS) -qarch=auto -qrealsize=8 -qdpc=e \
		-qspillsize=2500 -g -qfullpath $(QCMP)
FREEFLAGS  := -qsuffix=f=f90:cpp=F90
FIXEDFLAGS := -qfixed=132
#
# As of Apr/15/2003 cam2_0_2_dev18 xlfrte8.1.0.3 using -lmass causes the code to core-dump
# when using DEBUG compiler options.
#
# ... hence the following ugly expedient ...
ifeq ($(DEBUG),TRUE)
  LDFLAGS := $(QLDR)
else
  LDFLAGS := $(QLDR)
#  LDFLAGS += -lmass
endif

# Use the thread-safe compilers even when disabling openMP.
ifeq ($(SPMD),TRUE)
  FC := mpxlf90_r
else
  FC := xlf90_r
  FFLAGS += -WF,-DHIDE_MPI
endif

ifeq ($(SMP),TRUE)
  # THREADED_PTHREADS is used by the timing library
  CFLAGS += -DTHREADED_PTHREADS
endif

ifeq ($(DEBUG),TRUE)
#
# Bounds checking is unreliable on the IBM.
# Sometimes you can get it to go if you turn threading off (by deleting -qsmp=omp)
# Only turn float-trapping on for debug mode as it's a 20% performance hit.
#
  FFLAGS += -qinitauto=FF911299 -qflttrap=ov:zero:inv:en
  ifeq ($(SMP),TRUE)
    FFLAGS += -qsmp=omp:noopt
    LDFLAGS += -qsmp=omp:noopt
  endif
else

# Check for override of default Fortran compiler optimizations
  ifeq ($(F_OPTIMIZATION_OVERRIDE),$(null))
# Inline when not debugging
    FORTRAN_OPTIMIZATION := -O3 -qstrict -Q
  endif
  FFLAGS += $(FORTRAN_OPTIMIZATION)
  ifeq ($(SMP),TRUE)
    FFLAGS += -qsmp=omp
    LDFLAGS += -qsmp=omp
  endif
endif

.SUFFIXES:
.SUFFIXES: .F .F90 .c .o

.F.o:
	$(FC) -c $(FIXEDFLAGS) $(FFLAGS) $<

.F90.o:
	$(FC) -c $(FREEFLAGS) $(FFLAGS) $<

.c.o:
	cc -c  $(CFLAGS) $<


endif

#------------------------------------------------------------------------
# SGI
#------------------------------------------------------------------------

ifeq ($(UNAMES),IRIX64)

# default is either pure SPMD or pure SMP parallelism
ifeq ($(SMP),$(null))
  ifeq ($(SPMD),TRUE)
    SMP := FALSE
  else
    SMP := TRUE
  endif
endif

ESMF_ARCH   := IRIX64
FC          := f90
CPP         := /lib/cpp

CPPFLAGS    := -P
CFLAGS      := $(cpp_path) -64 -DIRIX64 -O2 $(CPPDEF)
mod_path    := -I$(ESMF_MOD)/$(ESMF_ARCH) -I$(MOD_NETCDF)
FFLAGS      := $(cpp_path) $(mod_path) -64 -r8 -i4 -c -cpp -extend_source -DIRIX64 \
		$(CPPDEF)
LDFLAGS     = -64

ifeq ($(SMP),TRUE)
  # THREADED_OMP is used by the timing library
  CFLAGS += -DTHREADED_OMP
# Only enable OpenMP, not SGI specific parallelism
# Note that -lmp should precede -lmpi when running hybrid MPI/OpenMP
  FFLAGS  += -mp -MP:dsm=OFF -MP:old_mp=OFF
  LDFLAGS += -mp -lmp
endif

ifeq ($(SPMD),TRUE)
  FFLAGS  += -I$(INC_MPI)
  LDFLAGS += -L$(LIB_MPI) -lmpi
else
  FFLAGS  += -DHIDE_MPI
endif

# WARNING: -mp and -g together has been known to cause wrong answers
ifeq ($(DEBUG),TRUE)
  FFLAGS  += -g -DEBUG:trap_uninitialized=ON -C
else
# Check for override of default Fortran compiler optimizations
  ifeq ($(F_OPTIMIZATION_OVERRIDE),$(null))
    FORTRAN_OPTIMIZATION := -O2
  endif
  FFLAGS += $(FORTRAN_OPTIMIZATION)
endif

.SUFFIXES:
.SUFFIXES: .F .F90 .c .o

.F.o:
	$(FC) $(FFLAGS) $<

.F90.o:
	$(FC) $(FFLAGS) $<
.c.o:
	cc -c  $(cpp_path) $(CFLAGS) $<


endif

#------------------------------------------------------------------------
# SUN
#------------------------------------------------------------------------

ifeq ($(UNAMES),SunOS)
ESMF_ARCH = solaris
FC	:= f90
FC77	:= f77
CPP     := /usr/ccs/lib/cpp
CFLAGS  := $(cpp_path) -DSUNOS $(CPPDEF)
mod_path:= -M$(ESMF_MOD)/$(ESMF_ARCH) -M$(MOD_NETCDF)
# Don't use OpenMP on Solaris as it currently causes problems
FFLAGS  := $(cpp_path) $(mod_path) -xs -stackvar -Qoption f90comp -r8const \
		-xtypemap=real:64,integer:32 -e -DSUNOS -DHIDE_MPI \
		$(CPPDEF)
LDFLAGS := -L/opt/SUNWspro/lib -lf77compat -openmp -fast
SPEC_FFLAGS := $(cpp_path) -r8 -i4 -c

ifeq ($(DEBUG),TRUE)
  FFLAGS += -g -dalign
  SPEC_FFLAGS += -g -dalign
else
# Check for override of default Fortran compiler optimizations
  ifeq ($(F_OPTIMIZATION_OVERRIDE),$(null))
# Inline code when not debugging -inline
# Use symbol table (-g) as make execution stable and can find out where problems are.
# Turn off aggressive optimization options (fsimple=1, fns=no) so the climate
# will be similar to climate on the IBM
    FORTRAN_OPTIMIZATION := -inline=%auto -fast -g -fsimple=1 -fns=no -O4
  endif
  FFLAGS += $(FORTRAN_OPTIMIZATION)
endif

ifeq ($(SPMD),TRUE)
  FFLAGS  += -I$(INC_MPI)
  LDFLAGS += -L$(LIB_MPI) -lmpich -lnsl -lsocket
else
  FFLAGS  += -DHIDE_MPI
endif

.SUFFIXES:
.SUFFIXES: .F .F90 .c .o

.F90.o:
	$(FC) -c $(FFLAGS) $<

endif

#------------------------------------------------------------------------
# Linux
#------------------------------------------------------------------------

ifeq ($(UNAMES),Linux)

ESMF_ARCH = linux

# if the compiler is not specified, choose PGI pgf90

ifeq ($(USER_FC),$(null))
  FC := pgf90
  FPGI := pgi
else
  FC := $(USER_FC)
endif

# USER_FC set to mpipgf90 is for LLNL MCR system. Define FPGI to correspond to
# pgf90 or mpipgf90. Define FLLNL to correspond to mpipgf90 or mpiifc.

ifeq ($(FPGI),$(null))
  ifeq ($(FC),mpipgf90)
    FPGI := pgi
    FLLNL := llnl
  else
    FPGI := notpgi
  endif
endif

# if USER_CC is set, use it.  Otherwise use pgcc if pgf90 is the Fortran compiler,
# and cc if not

ifeq ($(USER_CC),$(null))
  ifeq ($(FPGI),pgi)
    CC := pgcc
  else
    CC := cc
  endif
else
  CC := $(USER_CC)
endif

# Figure out the ESMF architecture.
ifeq ($(FPGI),pgi)
  ifeq ($(CC),pgcc)
    ESMF_ARCH = linux_pgi
  else
    ESMF_ARCH = linux_gnupgf90
  endif
endif

CFLAGS = $(cpp_path) -DLINUX $(CPPDEF)
ifeq ($(CC),pgcc)
  CFLAGS += -fast
endif

# pgf90
# -DPGF90 is for phcs and gauaw which normally use r16 arithmetic but is unavailable under pgf90

ifeq ($(FPGI),pgi)

  ifeq ($(SMP),$(null))
    SMP := TRUE
  endif
  mod_path  := -I$(ESMF_MOD)/$(ESMF_ARCH) -I$(MOD_NETCDF)
  FFLAGS    := $(cpp_path) $(mod_path) -r8 -i4 $(CPPDEF) -Mrecursive -Mdalign -Mextend -DLINUX \
	       -DPGF90
  ifeq ($FC),mpipgf90)
    $FFLAGS += -byteswapio
  endif
  FREEFLAGS := -Mfree
  LDFLAGS :=

  ifeq ($(DEBUG),TRUE)
    FFLAGS      += -g -Ktrap=fp -Mbounds
    SPEC_FFLAGS := $(FFLAGS) 
  else
    SPEC_FFLAGS := $(FFLAGS) 
# Check for override of default Fortran compiler optimizations
    ifeq ($(F_OPTIMIZATION_OVERRIDE),$(null))
#TBH:  override to -O0 due to failure of fv06 test
      FORTRAN_OPTIMIZATION := -fast
    endif
    FFLAGS += $(FORTRAN_OPTIMIZATION)
  endif

  ifeq ($(SMP),TRUE)
    FFLAGS += -mp
    LDFLAGS += -mp
  endif

endif

#--------------------------------------------------------------------------
# ifc - Intel compiler on Linux @ LLNL
#--------------------------------------------------------------------------
ifeq ($(FC),mpiifc)

  ifeq ($(SMP),$(null))
    SMP := TRUE
  endif
  mod_path  := -I$(ESMF_MOD)/$(ESMF_ARCH) -I$(MOD_NETCDF)
  FFLAGS    := $(cpp_path) $(mod_path) -extend_source -r8 -WB -tpp7 -Zp16 -ip $(CPPDEF) -DLINUX
  LDFLAGS      = -static -lelan -lelan3 -lrmscall -lelf -L/usr/local/intel/mkl/lib/32 -lmkl_def -lmkl_p4

  ifeq ($(DEBUG),TRUE)
    FFLAGS      += -g 
    SPEC_FFLAGS := $(FFLAGS) 
  else
    FFLAGS += O3
    SPEC_FFLAGS := $(FFLAGS) 
  endif

  ifeq ($(SMP),TRUE)
    FFLAGS += -mp
    LDFLAGS += -mp
  endif

endif


# lf95
#
# Note that as of lf95 version 6.1 threading does NOT work because of
# ridiculously small per thread stacksize limits.
#
#  -CcdRR8 is an undocumented flag which promotes only vars declared "real", not "real*4"
# --trace produces a call traceback on abort
# --trap causes code to stop on divide by zero or overflow exceptions
# --pca prevents overwriting constant arguments
# --staticlink enables transportability of the executable
# --chk for basic compiler checking (a,e,s,u,x)
# --chkglobal for global checking

ifeq ($(FC),lf95)
  ESMF_ARCH   := linux_lf95
  mod_path    := -I$(ESMF_MOD)/$(ESMF_ARCH) -I$(MOD_NETCDF)
  FFLAGS      := $(cpp_path) $(mod_path) -CcdRR8 $(CPPDEF) --trace --trap --wide -DLINUX
  SPEC_FFLAGS := $(FFLAGS)
  LDFLAGS     := --staticlink
  ifeq ($(DEBUG),TRUE)
#TBH:  this works    FFLAGS += -g --chk --pca
#TBH:  this FAILS    FFLAGS += -g --chk a,e,s,u,x --pca
    FFLAGS += -g --chk a,e,s,u --pca
  else
# Check for override of default Fortran compiler optimizations
    ifeq ($(F_OPTIMIZATION_OVERRIDE),$(null))
      FORTRAN_OPTIMIZATION := -O
    endif
    FFLAGS += $(FORTRAN_OPTIMIZATION)
  endif
endif

# FFC f90

ifeq ($(FC),f90)
  mod_path    := -I$(ESMF_MOD)/$(ESMF_ARCH) -I$(MOD_NETCDF)
  FFLAGS      := -CcdRR8 -Cpp -X9 -fw -Am -Wa,-W $(cpp_path) $(mod_path) $(CPPDEF) -DLINUX
  LDFLAGS      = 
  ifeq ($(DEBUG),TRUE)
    FFLAGS += -g -H aseu
  else
# Check for override of default Fortran compiler optimizations
    ifeq ($(F_OPTIMIZATION_OVERRIDE),$(null))
      FORTRAN_OPTIMIZATION := -Kfast,eval,fastlib,auto
    endif
    FFLAGS += $(FORTRAN_OPTIMIZATION)
  endif
endif

ifeq ($(FLLNL),$(null))
  LMPI := -lmpich
else
  LMPI := -lmpi
endif

# Flags common to all compilers

ifeq ($(SPMD),TRUE)
  FFLAGS  += -I$(INC_MPI)
  LDFLAGS += -L$(LIB_MPI) $(LMPI)
else
  FFLAGS  += -DHIDE_MPI
endif

ifeq ($(SMP),TRUE)
  # THREADED_OMP and THREADED_PTHREADS are used by the timing library
  ifeq ($(CC),pgcc)
    CFLAGS += -DTHREADED_OMP
  else
    CFLAGS += -DTHREADED_PTHREADS
  endif
endif

.SUFFIXES:
.SUFFIXES: .F .F90 .c .o

ifeq ($(FC),pgf90)
#
# To fix hanging problem when using sld dynamics, compile sgexx without "-fast"
#

sgexx.o: sgexx.F
	$(FC) -c $(SPEC_FFLAGS) $<
endif

ifeq ($(FC),lf95)
# lahey fails on binary_io due to writing wrap areas
binary_io.o: binary_io.F90
	$(FC) -c $(SPEC_FFLAGS) $<
wrap_nf.o: wrap_nf.F90
	$(FC) -c $(SPEC_FFLAGS) $<
wrap_mpi.o: wrap_mpi.F90
	$(FC) -c $(SPEC_FFLAGS) $<
endif

.F90.o:
	$(FC) -c $(FREEFLAGS) $(FFLAGS) $<
.F.o:
	$(FC) -c $(FFLAGS) $<
.c.o:
	$(CC) -c $(CFLAGS) $<

endif

#------------------------------------------------------------------------
# OSF1
#------------------------------------------------------------------------

ifeq ($(UNAMES),OSF1)

ESMF_ARCH     := alpha
CFLAGS        := $(cpp_path) -DOSF1 -O2 -omp $(CPPDEF)
FC            := f90
mod_path      := -I$(ESMF_MOD)/$(ESMF_ARCH) -I$(MOD_NETCDF)
FFLAGS        := $(cpp_path) $(mod_path) $(CPPDEF) -r8 -i4 -c -omp -automatic -fpe3 \
		 -check omp_bindings
FFLAGS_DOTF90 := -DOSF1 -free -fpe3
FFLAGS_DOTF   := -extend_source -omp -automatic
LDFLAGS       := -omp -lcxml

ifeq ($(SPMD),TRUE)
  FFLAGS  += -I$(INC_MPI) 
  LDFLAGS += -lmpi
else
  FFLAGS  += -DHIDE_MPI
endif

ifeq ($(SMP),TRUE)
  # THREADED_OMP is used by the timing library
  CFLAGS += -DTHREADED_OMP
endif

ifeq ($(DEBUG),TRUE)
  FFLAGS += -g3 -C
else
# Check for override of default Fortran compiler optimizations
  ifeq ($(F_OPTIMIZATION_OVERRIDE),$(null))
# Inline when not debugging
    FORTRAN_OPTIMIZATION := -O2 -inline speed
  endif
  FFLAGS += $(FORTRAN_OPTIMIZATION)
endif

ifeq ($(NO_SWITCH),$(null))
  NO_SWITCH := FALSE
endif
ifneq ($(NO_SWITCH),TRUE)
  LDFLAGS += -lelan
endif

.SUFFIXES:
.SUFFIXES: .F .F90 .c .o

.F.o:
	$(FC) $(FFLAGS) $(FFLAGS_DOTF) $<

.F90.o:
	$(FC) $(FFLAGS) $(FFLAGS_DOTF90) $<

.c.o:
	cc -c  $(CFLAGS) $<


endif

#------------------------------------------------------------------------
# Targets/rules that depend on architecture specific variables.
#------------------------------------------------------------------------

# The ESMF library is not made to be built in parallel, so specify only one job is to run.
$(ESMF_LIB)/$(ESMF_ARCH)/libesmf.a:
	cd $(ESMF_ROOT); \
           echo "Build the ESMF library."; \
           echo "ESMF is NOT supported by the CCSM project, but by the ESMF core team in NCAR/SCD"; \
           echo "See http://www.esmf.ucar.edu"; \
           $(MAKE) -j 1 BOPT=$(ESMF_BOPT) ESMF_BUILD=$(ESMF_BLD) ESMF_DIR=$(ESMF_ROOT) ESMF_ARCH=$(ESMF_ARCH);

time_manager.o : $(ESMF_LIB)/$(ESMF_ARCH)/libesmf.a

RM := rm
# Add user defined compiler and load flags if set, and replace FC if USER option set.
CFLAGS  += $(USER_CFLAGS)
FFLAGS  += $(USER_FFLAGS)
LDFLAGS  += $(USER_LDFLAGS)
ifneq ($(USER_FC),$(null))
FC := $(USER_FC)
endif

clean:
	$(RM) -rf esmf
	$(RM) -f Depends Srcfiles *.o *.mod *.stb *.f90 $(MODEL_EXEDIR)/$(EXENAME)

realclean: clean

include $(CURDIR)/Depends
