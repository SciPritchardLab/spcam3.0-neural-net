#
#       w2html.pm               Erik Kluzek             May/12/98
#                               NCAR/CGD/CSM
#                               erik@ncar.ucar.edu
#                               (303)497-1326
#
#	Subroutines to convert the MS-Word copy of the Users
#	Guide (output in HTML format) back to separate 
#      web-based HTML version.
#
#
use strict;
require "Web_File.pm";
package Word2HTML;

sub new {
#
#  Construct a new Word2HTML object
#
  my $proto = shift;
  my $file = shift;

  my $self = {};
  $self->{IN_FH}  = undef;            # Input file handle
  $self->{INFILE} = undef;            # Input filename
  $self->{WEB}    = undef;            # Web_File object
  bless( $self );
  return( $self );
}

#
# Methods
#

sub Web {
#
# Return the Web-Object contained in the Word2HTML object
#
  my $self = shift;

  my $web = $self->{WEB};
  return( $web );
}

sub IN_File {
#
# Return the input file-handle contained in the Word2HTML object
#
  my $self = shift;

  my $infile = $self->{INFILE};
  return( $infile );
}

sub IN_fh {
#
# Return the input file-handle contained in the Word2HTML object
#
  my $self = shift;

  my $in_fh = $self->{IN_FH};
  return( $in_fh );
}

sub open_in {
#
# Open the input file
#
  my $self = shift;
  my $infile = shift;
  my $file = shift;

  my $nm = ref($self). "::open_in";
  if ( ! defined($infile) || $infile eq "" ) {
    die "ERROR($nm) infile not given\n";
  }
  my $web = $self->{WEB};
  $self->{INFILE} = $infile;
  print "Realign comments and line lengths: \n";
  print "./perl_mod/realign.pl $infile\n";
  system( "./perl_mod/realign.pl $infile" );
  open( INFILE, "<$infile" ) || die "Could not open file: $infile";
  $self->{IN_FH} = \*INFILE;
  $web->restart( $file );
}

sub put_back_file_headers {
#
# Put back the file headers in each of the internal files
#
  my $self = shift;
 
  my $in_fh = $self->{IN_FH};
  my $web = $self->{WEB};
  while( <$in_fh> ) {
    if ( $web->find_begin_of_file( $_ ) ) {
      my $file = $web->file;
      open( FILE, ">$file" ) || die "Could not open file: $file";
      $web->begin_new_file( \*FILE );
      print FILE $_;
      while( <$in_fh> ) {
        #
        # Test for end of file, if not print current line to file
        #
        if ( $web->test_end_of_file( $_, $in_fh, \*FILE ) ) {
          last;
        }
        else {
          $web->find_headers( $_ );
          print FILE $_;
        }
      }
    }
  } 
  print "Hit end of input file: \n";
  close( $in_fh );
}



1   # to make use or require happy
