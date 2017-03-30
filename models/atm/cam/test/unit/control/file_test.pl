#!/usr/bin/env perl
#
# Simple unit test of the filenames module.
#
use 5.004;   # Use at least this version of perl
use Env qw(CASE DEBUG EXENAME MODEL_CFGDIR DYNAMICS MODEL_BLDDIR MODEL_EXEDIR
           MODEL_SRCDIR NAMELIST NEWBUILD PCNST PERGRO 
           PLEV PNATS RESOLUTION RUNTYPE SCRIPT_DIR SPMD);
use strict;
use Cwd;

# SCRIPT_DIR location of the main script
if ( defined($SCRIPT_DIR) ) {
  chdir( $SCRIPT_DIR );
} else {
  $SCRIPT_DIR = cwd( );
}
use lib "..", "../../../bld";   # List of where to look for the Perl modules

require "CAM_run_unit_test.pm";

my $unit = CAM_run_unit_test->new;
#
# unit-testing configuration settings
#
$unit->do_archivelog( "no" );  # Don't archive log files
$NEWBUILD = "TRUE";
$DEBUG = "TRUE";
$EXENAME = "fileunittest";
$CASE = "unittest";
my $ROOT_DIR  = "$SCRIPT_DIR";
$ROOT_DIR =~ s/\/test\/unit\/control//;
$MODEL_CFGDIR = "$ROOT_DIR/bld";   # Configure scripts directory
$MODEL_SRCDIR = "$ROOT_DIR/src";         # Model source directory

#
# Configuration settings
# (Most of the following don't really matter for the unit-test, but must be set)
#
$PERGRO = "FALSE";
$RUNTYPE = "initial";
$NAMELIST = "unittest.input";
$DYNAMICS = "eul";
$PCNST    =  1;
$PNATS    =  1;
$PLEV     =  26;
$RESOLUTION   = "T42";
#
# Build and run general test both with and without SPMD on
#
my @list;
foreach my $spmd ( ("TRUE", "FALSE") ) {
  #
  # Setup the directories, configure and build
  #
  $unit->setup_directories;
  $SPMD = $spmd;
  $unit->chdir( "$MODEL_BLDDIR" );
  # List of filenames to use in build
  @list = ( "filenames_test.F90", "filenames.F90",    "shr_kind_mod.F90", 
               "shr_sys_mod.F90",    "time_manager.F90", "endrun.F90", 
               "string_utils.F90",   "precision.F90",    "pmgrid.F90", 
               "dycore.F90",         "mpishorthand.F" );
  #
  # If SPMD on, then need wrap_mpi.F90 and the timing library files
  #
  if ( $SPMD =~ /TRUE/ ) {
    push( @list, ("wrap_mpi.F90", "spmdinit.F90", "spmd_phys.F90", 
                  "spmd_dyn.F90", "constituents.F90", "infnan.F90", "pspect.F90",
                  "comspe.F90", "ppgrid.F90", "physconst.F90", 
                  "shr_const_mod.F90" ) );
    my @timing_list = ( "f_wrappers.c", "get_cpustamp.c", "get_thread_num.c",
                        "t_error.c", "t_initialize.c", "t_pclstr.c",
                        "t_pr.c", "t_reset.c", "t_setoption.c", 
                        "t_stamp.c", "t_start.c", "t_stop.c" );
    push( @list, @timing_list );
  }
  $unit->configure( $SCRIPT_DIR, @list );
  $unit->machine_specs( "build" );
  $unit->make;
#
# Run general test
#
  $unit->chdir( "$MODEL_EXEDIR" );
  $unit->machine_specs( "run" );
  $unit->build_namelist( "general" );
  $unit->run( "Run the general test", "general.$SPMD.log" );
}
#
# Run fail tests
#
$unit->build_namelist( "list" );
my $file = "listtests";
$unit->run( "List the possible tests", "$file" );
my @tests = $unit->list_fail_tests( $file );
system( "/bin/rm failtest:*.log" );
foreach my $test ( @tests ) {
  $unit->build_namelist( "$test" );
  my $log = "$test.log";
  $unit->run( "Run the test: $test", "$test.log", "FAIL" );
}
print "Look at the log of the tests in: $MODEL_EXEDIR to ensure ok\n"
