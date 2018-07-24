#
#	CAM_run_unit_test.pm			Erik Kluzek
#
#	Perl5 Object module extending the CAM_run.pm module
#	to handle unit-testing of modules in the Community Atmospheric Model (CAM).
#
#	Methods:
#
#	new ---------------- Constructor
#	run ---------------- Extend run method to include option to set for a FAIL test.
#	exec --------------- Extend CAM_run exec method so that if FAIL is set
#				the method fails only if status is not a failure.
#	Srcfiles ----------- Create the Srcfiles file from the list of files entered.
#	configure ---------- Extend the configure method to create Srcfiles file, 
#                           from list entered.
#	build_namelist ----- Over-ride the build_namelist method to create a simple
#				input file with the information entered into the method.
#	list_fail_tests ---- List the fail tests from a unit-test.
#
#	$Id: CAM_run_unit_test.pm,v 1.1.2.1 2002/02/16 00:51:30 erik Exp $
#
use 5.004;   # Use at least this version of perl
use strict;
#use diagnostics;

package CAM_run_unit_test;
#
# Env variables to use the shorthand notation with
#
use Env qw(NAMELIST);
@CAM_run_unit_test::ISA = "CAM_run";
use CAM_run;
#
# Use all of the CAM_run methods and data
#

sub new {
#
# Constructor
#
  my $class = shift;

  my $self = $class->SUPER::new;
  $self->{'FAILTEST'} = "no";      # If this is a fail test or not
  $self->{'CHANGELOG'} = undef;    # Don't add ChangeLog to log-files
  bless( $self, $class );
  $self->do_clean( "yes" );    # Clean out old directories
  return( $self );
}

sub Srcfiles {
#
# Create the Srcfiles list of Source files
#
  my $self = shift;
  my @list = @_;

  open( SRCFILES, ">Srcfiles" ) || die "ERROR: Trouble opening Srcfiles file\n";
  foreach my $file ( @list ) {
     print SRCFILES "$file\n";
  }
  close( SRCFILES );
}

sub Filepath_append {
#
# Add the unit-test directory to the Filepath
#
  my $self = shift;
  my $dir = shift;

  open( FILEPATH, ">>Filepath" ) || die "ERROR:: Could not open Filepath\n";
  print FILEPATH "$dir\n";
  close( FILEPATH );
}

sub build_namelist {
#
# Over-ride the build_namelist method to create a simple input file with the
# information entered into the method.
#
  my $self = shift;
  my $input = shift;

  open( NL, ">$NAMELIST" ) || die "ERROR:: Could not open $NAMELIST\n";
  print NL "$input\n";
  close( NL );
}

sub configure {
#
# Run the configure script
#
  my $self = shift;
  my $unitdir = shift;
  my @list = @_;

  $self->SUPER::configure;
  $self->Srcfiles( @list );
  $self->Filepath_append( $unitdir );
}

sub run {
#
# Extend the run method to accept the FAIL option
#
  my $self = shift;
  my $desc = shift;
  my $logfile = shift;
  my $fail = shift;

  if ( defined($fail) || $fail eq "FAIL" ) {
    $self->{'FAILTEST'} = "yes";
  }
  $self->SUPER::run( $desc, $logfile );
  $self->{'FAILTEST'} = "no";
}

sub exec {
#
# Over-ride exec method so that if FAIL option set, will die if not
# an error rather than dying on an error
#
  my $self = shift;
  my $command = shift;
  my $die_msg = shift;
  my $echo = shift;

  if ( ! defined($command)  ) {
    die "Command not given to exec method";
  }
  print "$command\n";
  # Use system, so that results will be echoed out, if echo option set
  if ( defined($echo) && ($echo == "echo") ) {
    system( $command );
  # Otherwise, use backtics so that results of command will not be seen
  } else {
    `$command`;
  }
  #
  # Die on an error
  #
  if ( $? != 0 && (! $self->do( "FAILTEST" )) ) {
    if ( defined($die_msg) ) { print "$die_msg\n"; }
    $self->die( "ERROR:: Trouble executing: $command" );
    return;
  #
  # For fail tests only die if error did not happen
  #
  } elsif ( $? == 0 && $self->do( "FAIL") ) {
    if ( defined($die_msg) ) { print "$die_msg\n"; }
    $self->die( "ERROR:: A fail test did not terminate: $command" );
    return;
  }
}

sub list_fail_tests {
#
# List the tests from a unit-test
#
  my $self = shift;
  my $file = shift;

  open( LIST, "<$file" ) || die "ERROR:: Could not open list-file: $file";
  #
  # Loop until find the "general" test
  #
  while( defined($_ = <LIST>) ) {
    if ( $_ =~ /general/ ) { last; }
  }
  #
  # Now loop until end of file adding the rest of the tests
  #
  my @list;
  while( defined($_ = <LIST>) ) {
     if ( /^([a-zA-Z0-9:_-]+)/ ) {
       push( @list, $1 );
     }
  }
  close( LIST );
  if ( $#list < 0 ) { die "ERROR:: No tests listed in file: $file\n"; }
  return( @list );
}

1   # To make use or require happy
