#!/bin/csh
#=======================================================================
#
#  ncar_batch.csh
#
#  This is a batch submission script for test-model.pl for DAO
#  platforms and batch ques. 
#
#  Currently runs on tropic/dycore.gsfc.nasa.gov
#
# Note: make sure the number of nodes and shared-memory CPU's used given
# 	in the batch submission agrees with the settings below.
#
# Usage: as given below but either...
#
#	env SCRIPT_DIR=`pwd` qsub dao_batch.csh
#	env SCRIPT_DIR=`pwd` llsubmit dao_batch.csh
#
#-----------------------------------------------------------------------
# Batch options for machine with NQS batch system. (tropic, daley, NAS)
# submit with
#       env SCRIPT_DIR=`pwd` qsub dao_batch.csh
#-----------------------------------------------------------------------
#QSUB -q fesb          # Name of the que (CHANGE THIS if needed)
#QSUB -l mpp_p=8       # Maximum number of processes (CHANGE THIS if needed)
#QSUB -x               # Export all Environment variables
#QSUB -j testbatch.sgi.log  # Output job log to this file
#QSUB -eo              # Put standard error and standard out in same file
#QSUB                  # End of options
#-----------------------------------------------------------------------
# Batch options for machine with PBS batch system. (prospect)
# submit with
#       env SCRIPT_DIR=`pwd` qsub ncar_batch.csh
#-----------------------------------------------------------------------
#PBS -l ncpus=8
#PBS -l walltime=6:00:00
#PBS -l mem=6gb
#PBS -S /bin/csh
#PBS -V
#PBS -j eo
#PBS -N batch.pbs.log
# ------------------------------
# Export all Environment variables
#PBS -N batch.pbs.log
#-----------------------------------------------------------------------
# Batch options for machine with loadleveler (IBM SP) (blackforest)
# submit with
#        env SCRIPT_DIR=`pwd` llsubmit ncar_batch.csh
# On IBM with batch submission the number of nodes won't change
# even when the PE configuration is changed.
#-----------------------------------------------------------------------
# Name of the que (CHANGE THIS if needed)
# @ class       = csl_pr
# Number of nodes (CHANGE THIS if needed)
# @ node        = 8
# Switch to use (CHANGE THIS if needed)
# @ network.MPI = csss,not_shared,us
# @ output      = testbatch.aix.log
# @ error       = testbatch.aix.err
# @ node_usage  = not_shared
# @ job_type    = parallel
# @ tasks_per_node = 1
# Export all Environment variables
# @ environment = COPY_ALL
# @ queue
#
#=======================================================================
if ( ! $?SCRIPT_DIR )then
  echo "ERROR:: The SCRIPT_DIR env variable is not set\!";
  echo "   Set SCRIPT_DIR to the location of test-model.pl";
  echo "   This can be done as either:";
  echo "      env SCRIPT_DIR=`pwd` llsubmit $0";
  echo "      env SCRIPT_DIR=`pwd` qsub $0";
  exit;
endif
set OS = `uname -s`;
switch ( $OS )
  case AIX:
     setenv SPMD_NODES 8;
     echo "Set SPMD_NODES to $SPMD_NODES";
     breaksw;
  case IRIX64:
#
# tropic.gsfc.nasa.gov
# 
     setenv OMP_NUM_THREADS 4
     setenv TMPDIR /tmp
     breaksw;
  case Linux:
#
# dycore.gsfc.nasa.gov
# 
     setenv LIB_NETCDF /usr/local/netcdf/lib
     setenv INC_NETCDF /usr/local/netcdf/include
     setenv LIB_MPI    /usr/local/mpi/lib
     setenv INC_MPI    /usr/local/mpi/include
     setenv USER_CC    gcc
     breaksw;
  default:
    echo "Use default values for number of nodes and shared memory CPUs";
    exit;
endsw

setenv LANDMODEL CLM2
setenv CSMDATA /share/fvgcm/CAM/inputdata
setenv LAB "dao"
setenv DYNAMICS fv
cd $SCRIPT_DIR
./test-model.pl -skip fv:1-13
exit 0
