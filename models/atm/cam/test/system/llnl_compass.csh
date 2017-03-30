#!/bin/csh
#=======================================================================
#
#  llnl_compass.csh
#
#  This is a batch submission script for test-model.pl for LLNL's Compass
#  cluster of Compaq machines.  At other sites this file can be copied
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
#	env SCRIPT_DIR=`pwd` psub llnl_compass.csh
#
#-----------------------------------------------------------------------
# Batch options for machine with NQS batch system. (East)
# submit with
#       env SCRIPT_DIR=`pwd` psub llnl_compass.csh
#-----------------------------------------------------------------------
#PSUB -ln 4                 # Number of nodes to access (CHANGE IF NEEDED)
#PSUB -c compass            # Which machine to send it to (CHANGE IF NEEDED)
#PSUB -x                    # Export all Environment variables
#PSUB -o testbatch.osf.log  # Send log output to this file
#PSUB -noDFS                # DFS not needed
#PSUB -me                   # Notify me when done
#PSUB -tM 05:30             # Time limit
#PSUB -eo                   # Combine error and log output
#PSUB                       # End of options
#-----------------------------------------------------------------------
#=======================================================================
if ( ! $?SCRIPT_DIR )then
  echo "ERROR:: The SCRIPT_DIR env variable is not set\!";
  echo "   Set SCRIPT_DIR to the location of test-model.pl";
  echo "   This can be done as:";
  echo "      env SCRIPT_DIR=`pwd` psub $0";
  exit;
endif
set OS = `uname -s`;
switch ( $OS )
  case OSF1:
     setenv SPMD TRUE;
     setenv SPMD_NODES 4;
     setenv SHMEM_CPUS 4;
     echo "Set SPMD mode on and SPMD_NODES: $SPMD_NODES and SHMEM_CPUS: $SHMEM_CPUS";
     breaksw;
  default:
    echo "Use default values for number of nodes and shared memory CPUs";
    exit;
endsw

setenv LAB "llnl";
cd $SCRIPT_DIR;
./test-model.pl -nofail;
