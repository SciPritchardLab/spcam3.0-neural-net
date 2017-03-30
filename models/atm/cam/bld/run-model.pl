#!/usr/bin/env perl
#=======================================================================
#
#  This is a CAM stand-alone atm model test script
#
# Usage:
#
# perl run-model.pl
#
#-----------------------------------------------------------------------
# Batch options for machine with NQS
# submit with
#       env SCRIPT_DIR=`pwd` qsub run-model.pl
#-----------------------------------------------------------------------
#QSUB -l mpp_p=8
#QSUB -q fesb
#QSUB -x
#QSUB -eo
#QSUB
#-----------------------------------------------------------------------
# Batch options for machine with loadleveler (IBM SP)
# submit with
#       llsubmit run-model.pl
#-----------------------------------------------------------------------
# @ output = runmodel.aix.log
# @ error  = runmodel.aix.err
# @ node_usage = not_shared
# @ class = csl_pr
# @ job_type = parallel
# @ node = 2
# @ tasks_per_node = 16
# @ network.MPI = csss,not_shared,us
# @ queue
#
#=======================================================================
use strict;
use 5.004;   # Use at least this version of perl
#
# Use the shortened syntax for the following ENV variables
#
use Env qw(SCRIPT_DIR DEPGEN);
use Cwd;

# SCRIPT_DIR location of the main script
if ( defined($SCRIPT_DIR) ) {
  chdir( $SCRIPT_DIR );
} else {
  $SCRIPT_DIR = cwd( );
}
use lib ".";   # List of where to look for the Perl modules
require "CAM_run.pm";

my $cam = CAM_run->new;

$cam->setenv( "RESUB_YEAR", 1 );    # Resubmit until you reach this year
$cam->process_args;

#-----------------------------------------------------------------------
# Set directories
#  CAMROOT      = Location of the main root to CAM.
#-----------------------------------------------------------------------
  my $ROOT_DIR  = "$SCRIPT_DIR";
  $ROOT_DIR =~ s/\/models\/atm\/cam\/bld//;
  $cam->CFG->setcfg( "CAMROOT", "$ROOT_DIR" );          # CAM root
#  $cam->CFG->setcfg( "SMP", "FALSE" );
#-----------------------------------------------------------------------
# Configuration settings
#-----------------------------------------------------------------------
  if ( $cam->CFG->cfg( "DYNAMICS" ) eq "fv" ) {
    $cam->CFG->setcfg( "RESOLUTION", "2x2.5" ); 
  }
#
# Settings controlling CPU's and parallelization
# These are set in machine_specs. (uncomment if you want to control them here)
#
#  $cam->setenv( "SHMEM_CPUS",  2 );       # Number of shared-memory CPU's to use
#  $cam->setenv( "SPMD_NODES",  1 );       # Number of cluster nodes to use with MPI
#  $cam->setenv( "SPMD_CPUS_ON_NODE", 1 ); # number of on-node CPU's to use with MPI
#  $cam->CFG->setcfg( "SPMD",        "FALSE" ); # Turn SPMD mode on or off
#
# Namelist parameters. Namelist filename, run type caseid
# (By default do a 1 day initial, run followed by a restart run)
# (Then if do_resubmit set to "yes" continue to resubmit this simulation)
# (until you read year $RESUB_YEAR)
#
  $cam->do_resubmit( "no" );
  $cam->setenv( "CASE", "cambld");   # Case name
  my $RUNTYPE;
  if ( -f $cam->rest_pfile ) {     # Do restart run if restart file exists
    $RUNTYPE = "restart";
  } else {                         # Otherwise do an initial run
    $RUNTYPE = "initial";
  }
#
# To make code automatically resubmit remove the following two lines
# (Also change RESUB_YEAR above to the year you want to end on)
#
  $RUNTYPE = "initial";      # Always do an initial run
  $cam->do_resubmit( "no" );    # Don't resubmit after finished
#
# Namelist name
#
  $cam->setenv( "RUNTYPE", $RUNTYPE );
  my $NAMELIST = "atm.parm.$RUNTYPE"; # Filename of namelist
  $cam->setenv( "NAMELIST", $NAMELIST );
#
# ncdata_vers the version number to use for default ncdata files
# Use default value in DefaultCAMEXPNamelist.xml
#
  $cam->setenv( "NCDATA_VERS", undef );
# $cam->setenv( "NCDATA_VERS", 1 );   # Uncomment this line to explicitly set version
#-----------------------------------------------------------------------
#
# Setup the namelist
#
#-----------------------------------------------------------------------
#
#initialize associative arrays to handle namelist
#$main::CAMEXP{'key'} = value;
#to set non-standard namelist values (LSMEXP for the LSM namelist,
#and CLMEXP for the CLM namelist).
#
#If datasets not given default datasets will be used
#Resolution dependent parameters will also be set when namelist is built
#if not specified here. The values set here will override any defaults...
#
#String data such as filenames need to have be surrounded by backquoted
#single quotes (\'). For example:
#
# $main::CAMEXP{'ncdata'}  = "\'/CAM/csm/input/atm/cam1/SEP1.T42L26.112000.nc\'";
#
  %main::CAMEXP = {}; 
  %main::LSMEXP = {};
  %main::CLMEXP = {};
  $main::CAMEXP{'iyear_ad'}  = 1950;
  if ( $RUNTYPE eq "initial" ) {
    $main::CAMEXP{'nelapse'}   = -1;
  } else {
    $main::CAMEXP{'nelapse'}   = -30;
  }
  $main::CAMEXP{'inithist'} = "\'MONTHLY\'";

#-----------------------------------------------------------------------
#
# Typically you will not have to edit below this point!
#
#-----------------------------------------------------------------------
  $cam->setup_directories;         # Setup the build, run and source directories  (and softlinks)
#
# Configure/build model
#
  $cam->chdir( "MODEL_BLDDIR" );
  $cam->configure;                 # Run the configure script
  $cam->make;                      # Build the model executable
  $cam->chdir( "$SCRIPT_DIR" );    # Change directory back
  $cam->chdir( "MODEL_EXEDIR" );   # Change to exedir
  $cam->build_namelist;            # Build the model namelist
#
# Run the model configuration
#
  $cam->run_time_env;              # Setup the ENV variables needed to run
  $cam->run;                            # Run the model
#
# Resubmit
#
  $cam->setenv( "RUNTYPE", "restart" );
  $cam->chdir( "$SCRIPT_DIR" );         # Change to exedir
  $cam->resubmit;                       # Resubmit the model simulation if set
