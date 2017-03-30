#!/bin/csh
#=======================================================================
#
#  nersc_batch.csh
#
#  This is a batch submission script for test-model.pl for nersc
#  platforms and batch ques. At other sites this file can be copied
#  and customized for the batch ques available there. The things that most
#  likely will need to be changed are marked as "CHANGE THIS". Also the options
#  sent to "test-model.pl" below can be customized for how you would
#  like to run it.
#
# Note: make sure the number of nodes and shared-memory CPU's used given
# 	in the batch submission agrees with the settings below.
#
# Usage: as given below
#
#	env SCRIPT_DIR=`pwd` llsubmit nersc_batch.csh
#
#-----------------------------------------------------------------------
# Batch options for machine with loadleveler (IBM SP) (seaborg)
# submit with
#        env SCRIPT_DIR=`pwd` llsubmit nersc_batch.csh
# On IBM with batch submission the number of nodes won't change
# even when the PE configuration is changed.
#-----------------------------------------------------------------------
# Name of the que (CHANGE THIS if needed)
# @ class       = premium
# Number of nodes (CHANGE THIS if needed)
# @ node        = 8
# Switch to use (CHANGE THIS if needed)
# @ network.MPI = csss,shared,us
# @ shell = /usr/bin/csh
# @ output      = testbatch.aix.log
# @ error       = testbatch.aix.err
# @ node_usage  = not_shared
# @ job_type    = parallel
# @ tasks_per_node = 1
# Limit wall clock time (CHANGE THIS if needed)
# @ wall_clock_limit = 02:30:00 
# Notify via. email when job is done
# @ notification = complete 
# Export all Environment variables
# @ environment = COPY_ALL
# @ queue
#
#=======================================================================
if ( ! $?SCRIPT_DIR )then
  echo "ERROR:: The SCRIPT_DIR env variable is not set\!";
  echo "   Set SCRIPT_DIR to the location of test-model.pl";
  echo "   This can be done as:";
  echo "      env SCRIPT_DIR=`pwd` llsubmit $0";
  exit;
endif
set OS = `uname -s`;
switch ( $OS )
  case AIX:
     setenv SPMD_NODES 2;
     echo "Set SPMD_NODES to $SPMD_NODES";
     breaksw;
  default:
    echo "Use default values for number of nodes and shared memory CPUs";
    exit;
endsw

setenv LAB "nersc";
module load netcdf
echo $NETCDF_DIR
cd $SCRIPT_DIR;
./test-model.pl -nofail;
