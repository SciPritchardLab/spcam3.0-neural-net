static char *machineinfo = "    \n\
Libraries compiled on Tue May 14 13:35:37 EDT 2024 on br012.ib.bridges2.psc.edu   \n\
Machine characteristics: Linux br012.ib.bridges2.psc.edu 4.18.0-477.27.1.el8_8.x86_64 #1 SMP Thu Aug 31 10:29:22 EDT 2023 x86_64 x86_64 x86_64 GNU/Linux   \n\
-----------------------------------------  \n\
Using C compiler: mpicc   -O  -I/jet/home/jlin96/repositories/spcam3.0-neural-net/models/utils/esmf -I/jet/home/jlin96/repositories/spcam3.0-neural-net/models/utils/esmf/src/include -I/jet/home/jlin96/repositories/spcam3.0-neural-net/models/utils/esmf/build/linux_intel -I/jet/home/jlin96/repositories/spcam3.0-neural-net/models/utils/esmf/include   -I/jet/home/jlin96/repositories/spcam3.0-neural-net/models/utils/esmf/src/Infrastructure/mpiuni                          -D__SDIR__=''   \n\
C Compiler version:  \n\
gcc (GCC) 8.5.0 20210514 (Red Hat 8.5.0-18)  \n\
Copyright (C) 2018 Free Software Foundation, Inc.  \n\
This is free software; see the source for copying conditions.  There is NO  \n\
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  \n\
  \n\
C++ Compiler version:  \n\
gcc (GCC) 8.5.0 20210514 (Red Hat 8.5.0-18)  \n\
Copyright (C) 2018 Free Software Foundation, Inc.  \n\
This is free software; see the source for copying conditions.  There is NO  \n\
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  \n\
  \n\
Using Fortran compiler: mpiifort -O -I/jet/home/jlin96/repositories/spcam3.0-neural-net/models/utils/esmf -I/jet/home/jlin96/repositories/spcam3.0-neural-net/models/utils/esmf/src/include -I/jet/home/jlin96/repositories/spcam3.0-neural-net/models/utils/esmf/build/linux_intel -I/jet/home/jlin96/repositories/spcam3.0-neural-net/models/utils/esmf/include   -I/jet/home/jlin96/repositories/spcam3.0-neural-net/models/utils/esmf/src/Infrastructure/mpiuni                           \n\
Fortran Compiler version:  \n\
Error: Command line argument is needed!  \n\
Simple script to compile and/or link MPI programs.  \n\
Usage: mpiifort [options] <files>  \n\
----------------------------------------------------------------------------  \n\
The following options are supported:  \n\
   -fc=<name> | -f90=<name>  \n\
                   specify a FORTRAN compiler name: i.e. -fc=ifort  \n\
   -echo           print the scripts during their execution  \n\
   -show           show command lines without real calling  \n\
   -show_env       show environment variables  \n\
   -config=<name>  specify a configuration file: i.e. -config=ifort for mpif90-ifort.conf file  \n\
   -v              print version info of mpiifort and its native compiler  \n\
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
                   i.e -link_mpi=opt|opt_mt|dbg|dbg_mt  \n\
   -norpath        disable rpath for compiler wrapper of the Intel(R) MPI Library  \n\
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
   I_MPI_LINK      specify the version of the Intel(R) MPI Library  \n\
   I_MPI_DEBUG_INFO_STRIP  \n\
                   turn on/off the debug information stripping during static linking  \n\
----------------------------------------------------------------------------  \n\
-----------------------------------------  \n\
Using ESMF flags:                       \n\
-----------------------------------------  \n\
Using configuration flags:  \n\
-----------------------------------------  \n\
Using include paths: -I/jet/home/jlin96/repositories/spcam3.0-neural-net/models/utils/esmf -I/jet/home/jlin96/repositories/spcam3.0-neural-net/models/utils/esmf/src/include -I/jet/home/jlin96/repositories/spcam3.0-neural-net/models/utils/esmf/build/linux_intel -I/jet/home/jlin96/repositories/spcam3.0-neural-net/models/utils/esmf/include   -I/jet/home/jlin96/repositories/spcam3.0-neural-net/models/utils/esmf/src/Infrastructure/mpiuni    \n\
-----------------------------------------  \n\
Using ESMF directory: /jet/home/jlin96/repositories/spcam3.0-neural-net/models/utils/esmf  \n\
Using ESMF arch: linux_intel  \n\
------------------------------------------  \n\
Using C linker: mpicc  -O  -Wl,-rpath,/ocean/projects/atm200007p/jlin96/nnspreadtesting_good/specific/coupling_folder/cam_folder_debug/obj/esmf/lib/libO/linux_intel   \n\
Using Fortran linker: mpiifort -O -Wl,-rpath,/ocean/projects/atm200007p/jlin96/nnspreadtesting_good/specific/coupling_folder/cam_folder_debug/obj/esmf/lib/libO/linux_intel   \n\
Using libraries: -L/ocean/projects/atm200007p/jlin96/nnspreadtesting_good/specific/coupling_folder/cam_folder_debug/obj/esmf/lib/libO/linux_intel "; 
