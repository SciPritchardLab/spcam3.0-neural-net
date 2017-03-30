#! /bin/tcsh -f
#
#=======================================================================
#
#  run-pc.csh
#
#  Generic batch submission script for PC-linux using PBS.  
#
#-----------------------------------------------------------------------
# Batch options for machine with PBS batch system. (anchorage)
# Usage for Lahey compiler (default): 
#   qsub run-pc.csh
# Usage for pgf90 compiler with pgcc: 
#   env OPT_BLD=pgf90-pgcc qsub run-pc.csh
# Usage for pgf90 compiler with gcc: 
#   env OPT_BLD=pgf90-gcc qsub run-pc.csh
#-----------------------------------------------------------------------
#
# Name of the queue (CHANGE THIS if needed)
#PBS -q long
# Maximum number of processes (CHANGE THIS if needed)
#PBS -l nodes=4
# output file base name
#PBS -N run-pc
# Put standard error and standard out in same file
#PBS -j oe
# Export all Environment variables
#PBS -V
# End of options
#=======================================================================

set OS = `uname -s`;
switch ( $OS )
  case Linux:
     if ( ! $?PBS_JOBID ) then
       echo "${0}: ERROR::  This batch script must be submitted via PBS";
       echo "${0}:          on a Linux machine\!";
       exit;
     else
       echo "${0}: Running CAM on Linux using PBS";
     endif
     set job_id = `echo ${PBS_JOBID} | cut -f1 -d'.'`
     echo "${0}:  Set job_id to $job_id";
     if ( ! $?OPT_BLD ) then
       echo "${0}:  USING lf95 COMPILER"
       #TBH  Lahey stuff has HARD-CODED PATHS for now (gaaaak) ...
       echo "${0}:  which lf95 = `which lf95`"
       setenv USER_FC lf95
       echo "${0}:  Set USER_FC to $USER_FC";
       setenv INC_MPI /usr/local/mpich-lf95/include
       echo "${0}:  Set INC_MPI to $INC_MPI";
       setenv LIB_MPI /usr/local/mpich-lf95/lib
       echo "${0}:  Set LIB_MPI to $LIB_MPI";
       set mpirun = /usr/local/mpich-lf95/bin/mpirun
       echo "${0}:  Set mpirun to $mpirun";
       setenv INC_NETCDF /usr/local/netcdf-gcc-lf95/include
       echo "${0}:  Set INC_NETCDF to $INC_NETCDF";
       setenv MOD_NETCDF $INC_NETCDF
       echo "${0}:  Set MOD_NETCDF to $MOD_NETCDF";
       setenv LIB_NETCDF /usr/local/netcdf-gcc-lf95/lib
       echo "${0}:  Set LIB_NETCDF to $LIB_NETCDF";
     else
       if ( $OPT_BLD == "pgf90-gcc" ) then
         echo "${0}:  USING pgf90 COMPILER WITH gcc"
         setenv USER_CC gcc
         echo "${0}:  Set USER_CC to $USER_CC";
       else
         echo "${0}:  USING pgf90 COMPILER WITH pgcc"
       endif
       setenv INC_MPI /usr/local/mpich/include
       echo "${0}:  Set INC_MPI to $INC_MPI";
       setenv LIB_MPI /usr/local/mpich/lib
       echo "${0}:  Set LIB_MPI to $LIB_MPI";
       set mpirun = /usr/local/mpich/bin/mpirun
       echo "${0}:  Set mpirun to $mpirun";
       setenv INC_NETCDF /usr/local/netcdf/include
       echo "${0}:  Set INC_NETCDF to $INC_NETCDF";
       setenv MOD_NETCDF $INC_NETCDF
       echo "${0}:  Set MOD_NETCDF to $MOD_NETCDF";
       setenv LIB_NETCDF /usr/local/netcdf/lib
       echo "${0}:  Set LIB_NETCDF to $LIB_NETCDF";
     endif
     breaksw;
  default:
    echo "${0}:  Use default values for number of nodes and shared memory CPUs";    exit;
endsw

cd ${PBS_O_WORKDIR};

## Do our best to get sufficient stack memory
limit stacksize unlimited

## ROOT OF CAM DISTRIBUTION - probably needs to be customized.
## Contains the source code for the CAM distribution.
## (the root directory contains the subdirectory "models")
set camroot      = /fs/cgd/data0/$LOGNAME/cam2_0_2_devNN/cam1

## ROOT OF CAM DATA DISTRIBUTION - needs to be customized unless running at NCAR.
## Contains the initial and boundary data for the CAM distribution.
## (the root directory contains the subdirectories "atm" and "lnd")
setenv CSMDATA     /fs/cgd/csm/inputdata

## Default namelist settings:
## $case is the case identifier for this run. It will be placed in the namelist.
## $runtype is the run type: initial, restart, or branch.
## $nelapse is the number of timesteps to integrate, or number of days if negative.
set case         = camrun.$job_id
set runtype      = initial
set nelapse      = -1

## $wrkdir is a working directory where the model will be built and run.
## $blddir is the directory where model will be compiled.
## $rundir is the directory where the model will be run.
## $cfgdir is the directory containing the CAM configuration scripts.
set wrkdir       = /ptmp/$LOGNAME
set blddir       = $wrkdir/$case/bld
set rundir       = $wrkdir/$case
set cfgdir       = $camroot/models/atm/cam/bld

## Ensure that run and build directories exist
mkdir -p $rundir                || echo "cannot create $rundir" && exit 1
mkdir -p $blddir                || echo "cannot create $blddir" && exit 1

## If an executable doesn't exist, build one.
if ( ! -x $blddir/cam ) then
    cd $blddir                  || echo "cd $blddir failed" && exit 1
    $cfgdir/configure -spmd     || echo "configure failed" && exit 1
    echo "building CAM in $blddir ..."
    rm -f Depends
    gmake >&! MAKE.out      || echo "CAM build failed: see $blddir/MAKE.out" && exit 1
endif

## Create the namelist
cd $blddir                      || echo "cd $blddir failed" && exit 1
$cfgdir/build-namelist -s -case $case -runtype $runtype -o $rundir/namelist \
 -namelist "&camexp nelapse=$nelapse mss_irt=0 /"  || echo "build-namelist failed" && exit 1

## Run CAM
cd $rundir                      || echo "cd $rundir failed" && exit 1
echo "running CAM in $rundir"
$mpirun -np 4 $blddir/cam < namelist  || echo "CAM run failed" && exit 1

exit 0
