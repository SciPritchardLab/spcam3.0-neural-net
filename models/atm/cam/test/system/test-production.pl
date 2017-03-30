#!/usr/bin/env perl
#=======================================================================
#
#  This is a CAM stand-alone atm model test script
#  for doing a short several month test of the production
#  version of the model.
#
# Usage:
#
# perl test-production.pl [options] (dynamics)
#
#-----------------------------------------------------------------------
# Batch options for machine with loadleveler (IBM SP) (blackforest)
# submit with
#       llsubmit test-production.pl                     # for eul
#	env DYNAMICS="sld" llsubmit test-production.pl  # for sld
#	env DYNAMICS="fv"  llsubmit test-production.pl  # for fv
#-----------------------------------------------------------------------
# @ output = runprod.aix.log
# @ error  = runprod.aix.err
# @ node_usage = not_shared
# @ class = csl_pr
# @ job_type = parallel
# @ node = 8
# @ tasks_per_node = 1
# @ network.MPI = csss,not_shared,us
# @ environment = COPY_ALL
# @ ja_report = yes
# @ queue
#
#=======================================================================
use strict;
use 5.004;   # Use at least this version of perl
use Getopt::Long;
#
# Use the shortened syntax for the following ENV variables
#
use Env qw(SCRIPT_DIR DYNAMICS);
use Cwd;

sub get_files {
#
# Get the files needed to report the performance of a previous case
#
  my $cam = shift;
  my $filesub = shift;  # History filename to extract
  my $logfile = shift;  # Logfile

  my $nm = "main::get_files";
  #
  # Error checking
  #
  if ( !defined($filesub) || !defined($logfile) ) {
    die "ERROR($nm):: get_files needs the arguments filesub and logfile\n";
  }
  #
  # Do enough setup on the control case that we can get the needed files
  #
  my $case = $cam->env("CASE");
  $cam->{CASE} = $case;
  if ( ! defined($case) ) {
    die "ERROR($nm):: env variable CASE is not defined when get_files called\n";
  }
  #
  # Now get the relevent files
  #
  $cam->{'logfile'} = $logfile;
  my $pwd = cwd();
  my $MODEL_EXEDIR = $cam->CFG->cfg( "MODEL_EXEDIR" );
  $cam->chdir( $MODEL_EXEDIR );
  $cam->msrcp( "logs", $logfile, 1 );
  $cam->msrcp( "logs", "$MODEL_EXEDIR/timing.*", 1 );
  my $file = "$case.$filesub";
  $cam->msrcp( "hist", "$file", 1 );
  $cam->chdir( $pwd );
}

#---------------------------------------------------------------------------------------
# Main program
#---------------------------------------------------------------------------------------
#
# Parse the input arguments
#
  my $dyn = undef;
  GetOptions( "dyn=s" => \$dyn ) or die "ERROR:: Bad command line arguments to $0\n";
  if ( ! defined($dyn) && ($DYNAMICS =~ /./) ) {
    $dyn = $DYNAMICS;
  } else {
    $dyn = "eul";
  }
  print "Dy-core: $dyn\n";
#
# Information on control case
#
  my $CONTROL_CASE;
  my $CONTROL_NAME;
  if ( $dyn eq "eul" ) {
     $CONTROL_CASE = "eul_c2002-12-18_cam2_0_1_dev13";
     $CONTROL_NAME = "erik";
  } elsif ( $dyn eq "sld" ) {
     $CONTROL_CASE = undef;
  } elsif ( $dyn eq "fv" ) {
     $CONTROL_CASE = undef;
  }
#
# SCRIPT_DIR location of the main script
#
  if ( defined($SCRIPT_DIR) ) {
    chdir( $SCRIPT_DIR );
  } else {
    $SCRIPT_DIR = cwd( );
  }
  use lib ".", "../../bld";   # List of where to look for the Perl modules
  require "CAM_test.pm";
  require "cam_timing.pm";

  my $cam = CAM_test->new;

  $cam->CFG->setcfg( "DYNAMICS", $dyn );
  $cam->do_clean( "yes" );
#-----------------------------------------------------------------------
# Set directories
#  CAMROOT      = Location of the main root to CAM.
#-----------------------------------------------------------------------
  my $ROOT_DIR  = "$SCRIPT_DIR";
  $ROOT_DIR =~ s/\/models\/atm\/cam\/test\/system//;
  $cam->CFG->setcfg( "CAMROOT", "$ROOT_DIR" );          # CAM root

#-----------------------------------------------------------------------
# Configuration settings
#-----------------------------------------------------------------------
  if ( $cam->CFG->cfg( "DYNAMICS" ) eq "fv" ) {
    $cam->CFG->setcfg( "RESOLUTION", "2x2.5" ); 
  }
#
# Namelist parameters. Namelist filename, run type caseid
# (By default do a 1 day initial, run followed by a restart run)
#
  my $date = `date +"c%Y-%m-%d"`;
  chomp( $date );
  my $version = $cam->version;
  $cam->setenv( "CASE", "${dyn}_${date}_${version}");   # Case name
  my $RUNTYPE = "initial";      # Always do an initial run
  $cam->do_resubmit( "no" );    # Don't resubmit after finished
#
# Namelist name
#
  $cam->setenv( "RUNTYPE", $RUNTYPE );
  my $NAMELIST = "atm.parm.$RUNTYPE"; # Filename of namelist
  $cam->setenv( "NAMELIST", $NAMELIST );
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
  %main::CAMEXP  = {}; 
  %main::LSMEXP  = {};
  %main::CLMEXP  = {};
  %main::MPRUN2D = {};
  $main::CAMEXP{'ctitle'}    = "\'Standard production simulation for $dyn dynamics\'";
  $main::CAMEXP{'mss_irt'}   = 4096;
  $main::CAMEXP{'stop_ymd'}  = 10201;
  $main::CAMEXP{'fincl2'}    = "\'T:I\', \'PS:I\'";

  $cam->setup_directories;         # Setup the build, run and source directories  (and softlinks)
#
# Configure/build model
#
  $cam->chdir( "MODEL_BLDDIR" );
  $cam->configure;                 # Run the configure script
  $cam->make;                      # Build the model executable
  $cam->chdir( "MODEL_EXEDIR" );   # Change to exedir
  $cam->build_namelist( 1 );       # Build the model namelist
#
# Run the model configuration
#
  $cam->run_time_env;              # Setup the ENV variables needed to run
  my $log = "startup.log";
  my $log_file = "$log";
  $cam->run( "Startup simulation and run until:".$main::CAMEXP{'stop_ymd'}, $log_file );
  my $MODEL_EXEDIR = $cam->CFG->cfg( "MODEL_EXEDIR" );
#
# Store information needed for performance analysis
#
  foreach my $file ( glob("timing.*") ) {
    $cam->msrcp( "logs", "$file" );
  }
  my $CASE = $cam->env( "CASE" );
  my %opts = ( dir=>"$MODEL_EXEDIR", log=>"$log_file", case=>"$CASE" );
  my $cam_time = cam_timing->new( \%opts );
#
# Compare this simulation to previous control simulation
#
  my $MODEL_EXEDIR = $cam->CFG->cfg( "MODEL_EXEDIR" );
  my $filesub = "cam2.h0.0001-01.nc";
  my $file1 = "$MODEL_EXEDIR/$CASE.$filesub";
  $cam->msrcp( "hist", $file1, 1 );
  my $control = CAM_test->new;
  if ( defined($CONTROL_CASE) ) {
    print "Compare to previous control simulation\n";
    $control->setenv( "CASE", $CONTROL_CASE );
    $control->CFG->setcfg( "CAMROOT", "$ROOT_DIR" ); # set CAMROOT to something
    my $LOGNAME = $ENV{LOGNAME};
    $ENV{LOGNAME} = $CONTROL_NAME;
    $control->setup_directories;
    #
    # Get the config file from mss and read it in
    #
    $control->setenv( "RUNTYPE", $RUNTYPE );
    $control->setenv( "NAMELIST", $NAMELIST );
    my $MODEL_BLDDIR = $control->CFG->cfg( "MODEL_BLDDIR" );
    my $config_file = "$MODEL_BLDDIR/" . $control->config_file;
    $control->msrcp( "logs", $config_file, 1 );
    $control->CFG->read_config_cache( $config_file );
    #
    # Get relevent files from MSS
    #
    &get_files( $control, $filesub, "$log" );
    $ENV{LOGNAME} = $LOGNAME;
    my $MODEL_EXEDIR = $control->CFG->cfg( "MODEL_EXEDIR" );
    my $file2 = "$MODEL_EXEDIR/$CONTROL_CASE.$filesub";
    $cam->compare_files( $file1, $file2, "$MODEL_EXEDIR/control.cprout", 
                       "Not identical to control simulation: $CONTROL_CASE" );
  }
#
# Get information on the timing of the simulation versus the control simulation
#
  print "Performance analysis: \n";
  $cam_time->report_perf;
  if ( defined($CONTROL_CASE) ) {
    my $CASE = $control->env( "CASE" );
    my $LOG_DIR = $control->env( "LOG_DIR" );
    my $log_file = "$LOG_DIR/$log";
    my $MODEL_EXEDIR = $control->CFG->cfg( "MODEL_EXEDIR" );
    my %opts = ( dir=>"$MODEL_EXEDIR", log=>"$log_file", case=>"$CONTROL_CASE" );
    my $cont_time = cam_timing->new( \%opts );
    $cont_time->report_perf;
    $cam_time->compare_perf( $cont_time );
    $control->exec( "/bin/rm $MODEL_EXEDIR/*.nc *.log *.initial cam timing.*" );  # Delete files
  }
  $cam->chdir( "MODEL_EXEDIR" );
  $cam->exec( "/bin/rm *.nc *.log *.initial cam timing.*" );  # Delete files
