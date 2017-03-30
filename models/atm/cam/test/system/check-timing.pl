#!/usr/bin/env perl
#=======================================================================
#
#  This is a simple script to parse the timing files for a
#  CAM simulation in order to get gross timing information on
#  a simulation.
#
# Usage:
#
# perl check-timing.pl directory1 directory2
#
#=======================================================================
use English;
use strict;
use 5.004;   # Use at least this version of perl
use Cwd;
use Getopt::Long;
use lib (".", "../../bld");

#
# Get program name and directory path
#
my $ProgName;
($ProgName = $PROGRAM_NAME) =~ s!(.*)/!!;         # name of program
my $ProgDir = $1;                      # name of directory where program lives
# Add $ProgDir to the list of paths that Perl searches for modules
unshift @INC, "$ProgDir";
require "cam_timing.pm";

#
# Parse the input arguments
#
my %opts1 = ( log=>undef, dir=>undef, case=>undef, list=>undef );
my %opts2 = ( log=>undef, dir=>undef, case=>undef, list=>undef );
my $timer_list = undef;
my $help = undef;
GetOptions(
    "log1=s"    => \$opts1{log},
    "log2=s"    => \$opts2{log},
    "case1=s"   => \$opts1{case},
    "case2=s"   => \$opts2{case},
    "list=s"    => \$timer_list,
    "help"      => \$help
) or &usage( "Bad options" );
  if ( defined($help) ) {
     &usage( " " );
  }
  if ( defined($timer_list) ) {
    my @list = split(/,/, $timer_list );
    $opts1{list} = \@list;
    $opts2{list} = $opts1{list};
  }
  $opts1{dir} = shift(@main::ARGV);
  if ( ! defined($opts1{dir}) ) { &usage( "no arguments" ); }
  if ( ! -d "$opts1{dir}" ) { &usage( "directory1 not a directory: $opts1{dir}" ); }
  $opts2{dir} = shift(@main::ARGV);
  if ( defined($opts2{dir}) && (! -d $opts2{dir}) ) { &usage( "directory2 not a directory" ); }
  if (@ARGV) {
      &usage("unrecognized arguments: @ARGV" );
  }

#
# Get information on the timing of the whole simulation
#
  my $cam = cam_timing->new( \%opts1 );
  $cam->report_perf;
  if ( defined($opts2{dir}) ) {
    my $control = cam_timing->new( \%opts2 );
    $control->report_perf;
    print "\n\n";
    $cam->compare_perf( $control );
  }

sub usage {
#
# Report usage of script and die
#
  my $line = shift;

  if ( defined($line) ) {
     print "Usage error: $line\n";
  }
  die <<EOF;
SYNOPSIS
     $PROGRAM_NAME [options] directory1 directory2
OPTIONS
     directory1      Directory where CAM was run.
     directory2      (optional) another directory where CAM was run
                     (To compare the timing results of the first simulation to the second).
     -log1 <file>    (optional) Log-file for directory1.
     -log2 <file>    (optional) Log-file for directory2.
     -case1 <caseid> (optional) case name for directory1.
     -case2 <caseid> (optional) case name for directory2.
     -list <timers>  (optional) list of timers to extract times for.
     -help           List this help screen.

Giving log-file(s) allows the simulation length to be checked and to ensure
that the simulation(s) ran to completion.

example:

	$PROGRAM_NAME /ptmp/erik/case1 /ptmp/erik/case2 -list "dynpkg,physpkg"

	Will extract the timing information for the two cases as well as the
	timing information (maximum times) for dynpkg and physpkg.

EOF
}
