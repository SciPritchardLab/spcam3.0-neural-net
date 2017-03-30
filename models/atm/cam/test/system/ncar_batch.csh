#!/bin/csh
#=======================================================================
#
#  ncar_batch.csh
#
#  This is a batch submission script for test-model.pl for NCAR
#  platforms and batch queues. At other sites this file can be copied
#  and customized for the batch queues available there. The things that most
#  likely will need to be changed are marked as "CHANGE THIS". Also the options
#  sent to "test-model.pl" below can be customized for how you would
#  like to run it.
#
# Note: make sure the number of nodes and shared-memory CPU's used given
#       in the batch submission agrees with the settings below.
#
# Note: LOTS of refactoring is needed!!  There is WAY too much duplication.  
#
# Usage: as given below for each supported machine.
#
#-----------------------------------------------------------------------
# Batch options for machine with NQS batch system. (chinookfe)
# Usage: 
#   qsub ncar_batch.csh
# OR, to compare with a previous version whose source was checked out in
# directory /tmp/oldversion:  
#   env COMPARE_DIR=/tmp/oldversion/cam1 qsub ncar_batch.csh
#-----------------------------------------------------------------------
#QSUB -q ded_16        # Name of the queue (CHANGE THIS if needed)
#QSUB -l mpp_p=20      # Maximum number of processes (CHANGE THIS if needed)
#QSUB -x               # Export all Environment variables
#QSUB -eo              # Put standard error and standard out in same file
#QSUB -J y             # Put job log in its own file
#QSUB                  # End of options
#-----------------------------------------------------------------------
# Batch options for machine with PBS batch system. (anchorage)
# Usage for Lahey compiler (default): 
#   qsub ncar_batch.csh
# OR, to compare with a previous version whose source was checked out in
# directory /tmp/oldversion:  
#   env COMPARE_DIR=/tmp/oldversion/cam1 qsub ncar_batch.csh
# Usage for pgf90 compiler with pgcc: 
#   env OPT_BLD=pgf90-pgcc qsub ncar_batch.csh
# OR, to compare with a previous version whose source was checked out in
# directory /tmp/oldversion:  
#   env COMPARE_DIR=/tmp/oldversion/cam1 OPT_BLD=pgf90-pgcc qsub ncar_batch.csh
# Usage for pgf90 compiler with gcc: 
#   env OPT_BLD=pgf90-gcc qsub ncar_batch.csh
# OR, to compare with a previous version whose source was checked out in
# directory /tmp/oldversion:  
#   env COMPARE_DIR=/tmp/oldversion/cam1 OPT_BLD=pgf90-gcc qsub ncar_batch.csh
#-----------------------------------------------------------------------
# Name of the queue (CHANGE THIS if needed)
#PBS -q long
# Maximum number of processes (CHANGE THIS if needed)
#PBS -l nodes=4
# output file base name
#PBS -N testbatch.pbs
# Put standard error and standard out in same file
#PBS -j oe
# Export all Environment variables
#PBS -V
# End of options
#-----------------------------------------------------------------------
# Batch options for machine with loadleveler (IBM SP) (blackforest)
# Usage: 
#   llsubmit ncar_batch.csh
# OR, to compare with a previous version whose source was checked out in
# directory /tmp/oldversion:  
#   env COMPARE_DIR=/tmp/oldversion/cam1 llsubmit ncar_batch.csh
# On IBM with batch submission the number of nodes won't change
# even when the PE configuration is changed.
#-----------------------------------------------------------------------
# Name of the queue (CHANGE THIS if needed)
# bluesky
#  # @ class       = csl_rg8
# blackforest
# @ class       = csl_pr
# Number of nodes (CHANGE THIS if needed)
# @ node        = 8
# Switch to use (CHANGE THIS if needed)
# @ network.MPI = csss,not_shared,us
# @ output      = testbatch.aix.$(jobid).log
# @ error       = testbatch.aix.$(jobid).err
# @ node_usage  = not_shared
# @ job_type    = parallel
# @ tasks_per_node = 1
# Export all Environment variables
# @ environment = COPY_ALL
# @ queue
#
#=======================================================================
set OS = `uname -s`;
switch ( $OS )
  case AIX:
     if ( ! $?LOADL_JOB_NAME ) then
       echo "${0}: ERROR::  This batch script must be submitted via";
       echo "${0}:          LoadLeveler on an AIX machine\!";
       exit;
     else
       echo "${0}: Running test-model on AIX using LoadLeveler";
     endif
     if ( ! $?SCRIPT_DIR ) then
       setenv SCRIPT_DIR ${LOADL_STEP_INITDIR};
       echo "${0}:  Set SCRIPT_DIR to $SCRIPT_DIR";
     else
       echo "${0}:  Using SCRIPT_DIR = $SCRIPT_DIR";
     endif
     set job_id = `echo ${LOADL_JOB_NAME} | cut -f2 -d'.'`;
     echo "${0}:  Set job_id to $job_id";
     if ( ! $?CASE_DIR ) then
       set mycasebase = /ptmp/${LOGNAME}/test-model.${job_id}
     else
       set mycasebase = ${CASE_DIR}/test-model.${job_id}
       unsetenv CASE_DIR
     endif
     setenv SPMD_NODES 8;
     echo "${0}:  Set SPMD_NODES to $SPMD_NODES";
     # 64-bit netcdf libraries, override common .cshrc settings here...  
     setenv INC_NETCDF /usr/local/include
     echo "${0}:  Set INC_NETCDF to $INC_NETCDF";
     setenv LIB_NETCDF /usr/local/lib64/r4i4
     echo "${0}:  Set LIB_NETCDF to $LIB_NETCDF";
     unsetenv MP_PROCS;
     breaksw;
  case IRIX64:
     if ( ! $?QSUB_REQID ) then
       echo "${0}: ERROR::  This batch script must be submitted via NQS";
       echo "${0}:          on an IRIX machine\!";
       exit;
     else
       echo "${0}: Running test-model on IRIX using NQS";
     endif
     if ( ! $?SCRIPT_DIR ) then
       setenv SCRIPT_DIR ${QSUB_WORKDIR};
       echo "${0}:  Set SCRIPT_DIR to $SCRIPT_DIR";
     else
       echo "${0}:  Using SCRIPT_DIR = $SCRIPT_DIR";
     endif
     set job_id = `echo ${QSUB_REQID} | cut -f1 -d'.'`
     echo "${0}:  Set job_id to $job_id";
     if ( ! $?CASE_DIR ) then
       set mycasebase = /ptmp/${LOGNAME}/test-model.${job_id}
     else
       set mycasebase = ${CASE_DIR}/test-model.${job_id}
       unsetenv CASE_DIR
     endif
     setenv SHMEM_CPUS 8;
     echo "${0}:  Set SHMEM_CPUS to $SHMEM_CPUS";
     breaksw;
  case Linux:
     if ( ! $?PBS_JOBID ) then
       echo "${0}: ERROR::  This batch script must be submitted via PBS";
       echo "${0}:          on a Linux machine\!";
       exit;
     else
       echo "${0}: Running test-model on Linux using PBS";
     endif
     if ( ! $?SCRIPT_DIR ) then
       setenv SCRIPT_DIR ${PBS_O_WORKDIR};
       echo "${0}:  Set SCRIPT_DIR to $SCRIPT_DIR";
     else
       echo "${0}:  Using SCRIPT_DIR = $SCRIPT_DIR";
     endif
     set job_id = `echo ${PBS_JOBID} | cut -f1 -d'.'`
     echo "${0}:  Set job_id to $job_id";
     if ( ! $?CASE_DIR ) then
       set mycasebase = /scratch/cluster/${LOGNAME}/test-model.${job_id}
     else
       set mycasebase = ${CASE_DIR}/test-model.${job_id}
       unsetenv CASE_DIR
     endif
     setenv SMP FALSE;
     echo "${0}:  Set SMP to $SMP";
     setenv SPMD TRUE;
     echo "${0}:  Set SPMD to $SPMD";
     setenv SPMD_NODES 4;
     echo "${0}:  Set SPMD_NODES to $SPMD_NODES";
     # TBH:  note that the value of $OPT_BLD is NOT checked yet...
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
       setenv SPMD_RUNCMND "/usr/local/mpich-lf95/bin/mpirun -np $SPMD_NODES"
       echo "${0}:  Set SPMD_RUNCMND to $SPMD_RUNCMND";
       setenv INC_NETCDF /usr/local/netcdf-gcc-lf95/include
       echo "${0}:  Set INC_NETCDF to $INC_NETCDF";
       setenv MOD_NETCDF $INC_NETCDF
       echo "${0}:  Set MOD_NETCDF to $MOD_NETCDF";
       setenv LIB_NETCDF /usr/local/netcdf-gcc-lf95/lib
       echo "${0}:  Set LIB_NETCDF to $LIB_NETCDF";
       set mycasebase = ${mycasebase}.lf95
     else
       if ( $OPT_BLD == "pgf90-gcc" ) then
         echo "${0}:  USING pgf90 COMPILER WITH gcc"
         setenv USER_CC gcc
         echo "${0}:  Set USER_CC to $USER_CC";
         set mycasebase = ${mycasebase}.pgf90-gcc
       else
         echo "${0}:  USING pgf90 COMPILER WITH pgcc"
         set mycasebase = ${mycasebase}.pgf90-pgcc
       endif
       setenv INC_MPI /usr/local/mpich/include
       echo "${0}:  Set INC_MPI to $INC_MPI";
       setenv LIB_MPI /usr/local/mpich/lib
       echo "${0}:  Set LIB_MPI to $LIB_MPI";
       setenv INC_NETCDF /usr/local/netcdf/include
       echo "${0}:  Set INC_NETCDF to $INC_NETCDF";
       setenv MOD_NETCDF $INC_NETCDF
       echo "${0}:  Set MOD_NETCDF to $MOD_NETCDF";
       setenv LIB_NETCDF /usr/local/netcdf/lib
       echo "${0}:  Set LIB_NETCDF to $LIB_NETCDF";
     endif
     breaksw;
  case OSF1:
     echo "${0}:  ERROR:  This script does not support OSF1 yet on NCAR";
     echo "${0}:          machines";
     exit;
     setenv SPMD TRUE;
     echo "${0}:  Set SPMD to $SPMD";
     setenv SPMD_NODES 4;
     echo "${0}:  Set SPMD_NODES to $SPMD_NODES";
     setenv SHMEM_CPUS 4;
     echo "${0}:  Set SHMEM_CPUS to $SHMEM_CPUS";
     breaksw;
  default:
    echo "${0}:  Use default values for number of nodes and shared memory CPUs";
    exit;
endsw
setenv CASE_DIR ${mycasebase}
echo "${0}:  Set CASE_DIR to $CASE_DIR";

setenv LAB "ncar";
cd $SCRIPT_DIR;

if ( ! $?COMPARE_DIR ) then
  # default:  no comparison with a previous version
  # NOTE:  -nofail removed to simplify debugging
  # ./test-model.pl -unique_id $job_id -nofail;
  ./test-model.pl -unique_id $job_id;
else
  # compare with a previous version
  echo "${0}:  NOTE:: The COMPARE_DIR env variable is set to $COMPARE_DIR";
  echo "${0}:  running test-model.pl with option < -compare $COMPARE_DIR >";
  # NOTE:  -nofail removed to simplify debugging
  # ./test-model.pl -unique_id $job_id -compare $COMPARE_DIR -nofail;
  ./test-model.pl -unique_id $job_id -compare $COMPARE_DIR;
endif

