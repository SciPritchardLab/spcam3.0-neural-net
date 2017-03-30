#!/usr/bin/env perl
#=======================================================================
#
#  This script will check that the namelist you have created will work.
#  It checks that the format is correct, that there are no duplicate
#  entries, and that needed setttings are given. If it can successfully
#  read the namelist it will print out a suggested one (which may have
#  some changes to it).
#
# Usage:
#
# perl check-namelist.pl filename
#
#=======================================================================
#
# Use the shortened syntax for the following ENV variables
#
use Env qw(AM_DATDIR CASE CASE_DIR DYNAMICS 
           LAB LM_DATDIR LOGNAME 
           MODEL_CFGDIR MODEL_DATDIR 
           MODEL_EXEDIR MODEL_SRCDIR NAMELIST PERGRO 
           RESOLUTION RUNTYPE 
           SCRIPT_DIR);
use strict;
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

$DYNAMICS = "eul"; # Model dynamics (over-ridden if commandline option given)
my $filename = shift( @ARGV );
if ( ! defined( $filename ) ) {
  die "ERROR:: Did not send filename: usage $0 filename\n";
}
$cam->config_env;  # Set the env variables from the last configure
$cam->do_interactive( "yes" );   # Turn interactive mode on
$cam->machine_specs( "run" );    # Setup the ENV variables specific to each machine-typ

$cam->Namelist->parse( $filename );    # Parse the namelist
