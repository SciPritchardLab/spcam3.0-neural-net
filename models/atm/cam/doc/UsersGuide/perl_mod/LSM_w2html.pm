#
#       LSM_w2html.pm           Erik Kluzek             May/18/98
#                               NCAR/CGD/CSM
#                               erik@ncar.ucar.edu
#                               (303)497-1326
#
#	Subroutines to convert the MS-Word copy of the LSM as
#	coupled to CCM Users Guide (output in HTML format) back to 
#	separate web-based HTML version.
#
#
use strict;
package LSM_Word2HTML;
@LSM_Word2HTML::ISA = "Word2HTML";
require "w2html.pm";
use LSM_Web_File;

#-----------------------------------------------------------
# Constructor
#-----------------------------------------------------------

sub new {
#
#  Construct a new LSM_Word2HTML object
#
  my $proto = shift;
  my $file = shift;

  my $self = $proto->SUPER::new( );
  $self->{WEB}    = LSM_Web_File->new( ); # Web_File object
  bless( $self );
  return( $self );
}

1   # to make use or require happy
