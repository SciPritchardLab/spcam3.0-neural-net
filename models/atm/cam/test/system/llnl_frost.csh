#!/bin/csh
#=======================================================================
#
#  llnl_frost.csh
#
#  This is a batch submission script for test-model.pl for LLNL's Frost Pacific
#  batch ques. At other sites this file can be copied
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
#	env SCRIPT_DIR=`pwd` psub llnl_frost.csh
#
#-----------------------------------------------------------------------
# Batch options for machine with NQS batch system. (Frost)
# submit with
#       env SCRIPT_DIR=`pwd` psub llnl_frost.csh
#-----------------------------------------------------------------------
#PSUB -c pbatch             # Batch pool to use (CHANGE IF NEEDED)
#PSUB -b climate            # Bank to use (CHANGE IF NEEDED)
#PSUB -ln 1 -g 4                 # Number of nodes to access (CHANGE IF NEEDED)
#PSUB -s /bin/csh           # Shell to use
#PSUB -x                    # Export all Environment variables
#PSUB -o testbatch.frost.log # Send log output to this file
#PSUB -c frost              # Constrain to frost machine
#PSUB -me                   # Notify me when done
#PSUB -tM 08:00             # Time limit
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
  case AIX:
     setenv SPMD_NODES 4;
     echo "Set SPMD_NODES to $SPMD_NODES";
     setenv SHMEM_CPUS 4
     echo "Set SHMEM_CPUS to $SHMEM_CPUS";
     setenv SPMD TRUE
     setenv CSMDATA /p/gf1/mirin/cam2data/inputdata
     breaksw;
  default:
    echo "Use default values for number of nodes and shared memory CPUs";
    exit;
endsw

setenv LAB "llnl";
cd $SCRIPT_DIR;
./test-model.pl -nofail;
