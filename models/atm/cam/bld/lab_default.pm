#
#	lab_default.pm			Erik Kluzek
#
#	Object to handle default settings where default
#	values are dependent not only on the lab but also
#	on the machine platform.
#
#	Methods:
#
#	new ----------------------- Constructor.
#	usage --------------------- If new isn't invoked correctly.
#	value --------------------- Return value for given Lab and O/S.
#
#	Sequential functions:
#
#	expand_env_vars_in_string - Expand any environment variables in the string.
#
#	$Id: lab_default.pm,v 1.4.10.2 2002/05/16 21:08:10 erik Exp $
#
use strict;
#use diagnostics;
package lab_default;

sub new {
#
# Constructor
#
  my $class = shift;
  my $lab_array_ref = shift;
  my @array = @_;
  my $self = {};

  my $set = undef;
  my $elem;
  #
  # Get lab list and convert it to a hash
  #
  my @labs = @$lab_array_ref;
  my %labs;
  foreach $elem ( @labs ) {
    $labs{$elem} = 1;
  }
  $self->{'LABLIST'} = $lab_array_ref;
  while ( defined($elem = shift( @array )) ) {
    if ( !defined($labs{$elem}) ) { &lab_default::usage( $class, "Lab $elem not in the lab list @labs" ); }
    my $hash = shift( @array );
    if ( ! defined($hash) ) { &lab_default::usage( $class, "Lab $elem not followed by anything" ); }
    if ( ref($hash) ne "HASH" ) { &lab_default::usage( $class, "$elem not followed by a hash reference" ); }
    my %hash = %$hash;
    if ( ! defined($hash{'default'} ) ) { &lab_default::usage( $class, "default not set for lab $elem" ); }
    $self->{$elem} = $hash;
    if ( $elem eq "default" ) { $set = 1; }
  }
  if ( ! defined($set) ) { &lab_default::usage( $class, "default lab not set" ); }
  bless( $self, $class );
  return( $self );
}

sub usage {
#
# Print out usage if something went wrong on the constructor call
#
  my $class = shift;
  my $prompt = shift;

  print "$prompt\n\n";
  print "Usage: \$value = $class->new( (\'LAB\', \\\%hash, \'LAB\', \\\%hash, ...) ) \n";
  print "    Error in the constructor call of $class\n";
  print "    Must be a series of different strings to describe a lab, followed by a hash\n";
  print "    Also \"default\" must be one of the lab descriptors\n";
  print "    Each hash should contain a default value as well as values for various platforms\n";
  die "Error on $class constructor";
}

sub expand_env_vars_in_string {
#
# Expand any environment variables that are in a string
# C-shell definition of variables is: 
#
# 	"A variable name consists of up to 20 letters and digits, 
# 	and starts with a  letter  (the underscore is considered 
#	a letter)."
#
  my $CAM = shift;
  my $value = shift;

  while ( $value =~ /^(.*)\${*([a-zA-Z_]+[a-zA-Z0-9_]{0,19})}*(.*)$/ ) {
    my $env = $2;
    my $lead = $1;
    my $tail = $3;
    my $env_value;
    # If CAM object has 
    if ( defined( $CAM ) ) {
      if ( $CAM->exists( $env ) ) {
        $env_value = $CAM->env($env);
        if ( ! defined( $env_value ) ) {
          die "ERROR:: Environment variable $env needed in setting this value $value";
        }
      }
      if ( $CAM->CFG->exists( $env ) ) {
        $env_value = $CAM->CFG->cfg($env);
        if ( ! defined( $env_value ) ) {
          die "ERROR:: Configuration variable $env needed in setting this value $value";
        }
      }
    } else {
      die "ERROR:: Environment or config variable $env needed in setting this value $value";
    }
    $value = "${lead}${env_value}${tail}";
  }
  return( $value );
}

sub value {
#
# Get the value for the given LAB and platform
#
  my $self = shift;
  my $lab = shift;
  my $OS = shift;
  my $CAM = shift;  # Optional CAM object with environment variables to use

  my $value = undef;
  if ( ! defined($self->{$lab} ) ) { $lab = "default"; }
  my $ref = $self->{$lab}; my %ref = %$ref;
  if ( ! defined($ref{$OS}) ) {
    $value = $ref{'default'};
  } else {
    $value = $ref{$OS};
  }
  # If it contains environment variables within it
  $value = &lab_default::expand_env_vars_in_string( $CAM, $value );
  return( $value );
}

1 # to make use or require happy
