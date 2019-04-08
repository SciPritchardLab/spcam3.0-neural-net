static char *machineinfo = "    \n\
Libraries compiled on Mon Apr  8 05:54:19 CDT 2019 on login1.stampede2.tacc.utexas.edu   \n\
Machine characteristics: Linux login1.stampede2.tacc.utexas.edu 3.10.0-957.5.1.el7.x86_64 #1 SMP Fri Feb 1 14:54:57 UTC 2019 x86_64 x86_64 x86_64 GNU/Linux   \n\
-----------------------------------------  \n\
Using C compiler: mpicc  -openmp  -O  -I/home1/05488/tg847872/repositories/spcam3.0-neural-net/models/utils/esmf -I/home1/05488/tg847872/repositories/spcam3.0-neural-net/models/utils/esmf/src/include -I/home1/05488/tg847872/repositories/spcam3.0-neural-net/models/utils/esmf/build/linux_intel -I/home1/05488/tg847872/repositories/spcam3.0-neural-net/models/utils/esmf/include   -I/home1/05488/tg847872/repositories/spcam3.0-neural-net/models/utils/esmf/src/Infrastructure/mpiuni                          -D__SDIR__=''   \n\
C Compiler version:  \n\
icc: command line remark #10411: option '-openmp' is deprecated and will be removed in a future release. Please use the replacement option '-qopenmp'  \n\
icc (ICC) 17.0.4 20170411  \n\
Copyright (C) 1985-2017 Intel Corporation.  All rights reserved.  \n\
  \n\
C++ Compiler version:  \n\
icc (ICC) 17.0.4 20170411  \n\
Copyright (C) 1985-2017 Intel Corporation.  All rights reserved.  \n\
  \n\
Using Fortran compiler: mpif90 -openmp  -O -I/home1/05488/tg847872/repositories/spcam3.0-neural-net/models/utils/esmf -I/home1/05488/tg847872/repositories/spcam3.0-neural-net/models/utils/esmf/src/include -I/home1/05488/tg847872/repositories/spcam3.0-neural-net/models/utils/esmf/build/linux_intel -I/home1/05488/tg847872/repositories/spcam3.0-neural-net/models/utils/esmf/include   -I/home1/05488/tg847872/repositories/spcam3.0-neural-net/models/utils/esmf/src/Infrastructure/mpiuni                           \n\
Fortran Compiler version:  \n\
Error: Command line argument is needed!  \n\
Simple script to compile and/or link MPI programs.  \n\
Usage: mpif90 [options] <files>  \n\
----------------------------------------------------------------------------  \n\
The following options are supported:  \n\
   -fc=<name> | -f90=<name>  \n\
                   specify a FORTRAN compiler name: i.e. -fc=ifort  \n\
   -echo           print the scripts during their execution  \n\
   -show           show command lines without real calling  \n\
   -config=<name>  specify a configuration file: i.e. -config=ifort for mpif90-ifort.conf file  \n\
   -v              print version info of mpif90 and its native compiler  \n\
   -profile=<name> specify a profile configuration file (an MPI profiling  \n\
                   library): i.e. -profile=myprofile for the myprofile.cfg file.  \n\
                   As a special case, lib<name>.so or lib<name>.a may be used  \n\
                   if the library is found  \n\
   -check_mpi      link against the Intel(R) Trace Collector (-profile=vtmc).  \n\
   -static_mpi     link the Intel(R) MPI Library statically  \n\
   -mt_mpi         link the thread safe version of the Intel(R) MPI Library  \n\
   -ilp64          link the ILP64 support of the Intel(R) MPI Library  \n\
   -no_ilp64       disable ILP64 support explicitly  \n\
   -fast           the same as -static_mpi + pass -fast option to a compiler.  \n\
   -t or -trace  \n\
                   link against the Intel(R) Trace Collector  \n\
   -trace-imbalance  \n\
                   link against the Intel(R) Trace Collector imbalance library  \n\
                   (-profile=vtim)  \n\
   -dynamic_log    link against the Intel(R) Trace Collector dynamically  \n\
   -static         use static linkage method  \n\
   -nostrip        turn off the debug information stripping during static linking  \n\
   -O              enable optimization  \n\
   -link_mpi=<name>  \n\
                   link against the specified version of the Intel(R) MPI Library  \n\
All other options will be passed to the compiler without changing.  \n\
----------------------------------------------------------------------------  \n\
The following environment variables are used:  \n\
   I_MPI_ROOT      the Intel(R) MPI Library installation directory path  \n\
   I_MPI_F90 or MPICH_F90  \n\
                   the path/name of the underlying compiler to be used  \n\
   I_MPI_FC_PROFILE or I_MPI_F90_PROFILE or MPIF90_PROFILE  \n\
                   the name of profile file (without extension)  \n\
   I_MPI_COMPILER_CONFIG_DIR  \n\
                   the folder which contains configuration files *.conf  \n\
   I_MPI_TRACE_PROFILE  \n\
                   specify a default profile for the -trace option  \n\
   I_MPI_CHECK_PROFILE  \n\
                   specify a default profile for the -check_mpi option  \n\
   I_MPI_CHECK_COMPILER  \n\
                   enable compiler setup checks  \n\
   I_MPI_LINK      specify the version of the Intel(R) MPI Library  \n\
   I_MPI_DEBUG_INFO_STRIP  \n\
                   turn on/off the debug information stripping during static linking  \n\
----------------------------------------------------------------------------  \n\
-----------------------------------------  \n\
Using ESMF flags:                       \n\
-----------------------------------------  \n\
Using configuration flags:  \n\
-----------------------------------------  \n\
Using include paths: -I/home1/05488/tg847872/repositories/spcam3.0-neural-net/models/utils/esmf -I/home1/05488/tg847872/repositories/spcam3.0-neural-net/models/utils/esmf/src/include -I/home1/05488/tg847872/repositories/spcam3.0-neural-net/models/utils/esmf/build/linux_intel -I/home1/05488/tg847872/repositories/spcam3.0-neural-net/models/utils/esmf/include   -I/home1/05488/tg847872/repositories/spcam3.0-neural-net/models/utils/esmf/src/Infrastructure/mpiuni    \n\
-----------------------------------------  \n\
Using ESMF directory: /home1/05488/tg847872/repositories/spcam3.0-neural-net/models/utils/esmf  \n\
Using ESMF arch: linux_intel  \n\
------------------------------------------  \n\
Using C linker: mpicc  -O  -Wl,-rpath,/home1/05488/tg847872/spcam3-nn-build/obj/esmf/lib/libO/linux_intel   \n\
Using Fortran linker: mpif90  -O -Wl,-rpath,/home1/05488/tg847872/spcam3-nn-build/obj/esmf/lib/libO/linux_intel   \n\
Using libraries: -L/home1/05488/tg847872/spcam3-nn-build/obj/esmf/lib/libO/linux_intel "; 
