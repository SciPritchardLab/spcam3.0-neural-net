#
#       CSM_w2html.pm           Erik Kluzek             May/18/98
#                               NCAR/CGD/CSM
#                               erik@ncar.ucar.edu
#                               (303)497-1326
#
#	Subroutines to convert the MS-Word copy of the CSM 
#	Users Guide (output in HTML format) back to 
#	separate web-based HTML version.
#
#
use strict;
package CSM_Word2HTML;
@CSM_Word2HTML::ISA = "Word2HTML";
require "w2html.pm";
use CSM_Web_File;

#-----------------------------------------------------------
# Constructor
#-----------------------------------------------------------

sub new {
#
#  Construct a new CSM_Word2HTML object
#
  my $proto = shift;
  my $file = shift;

  my $self = {};
  my $self = $proto->SUPER::new( );
  $self->{WEB}    = CSM_Web_File->new( ); # Web_File object
  bless( $self );
  return( $self );
}

1   # to make use or require happy
