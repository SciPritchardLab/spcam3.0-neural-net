#
#       html2w.pm               Erik Kluzek             Mar/2/98
#                               NCAR/CGD/CSM
#                               erik@ncar.ucar.edu
#                               (303)497-1326
#
#	Subroutines to convert the CCM users guide HTML source files to a file 
#	that can be read in Word for a paper-based copy.
#
#
use strict;
require "Web_File.pm";
package HTML2Word;

sub new {
#
#  Construct a new HTML2Word object
#
  my $proto = shift;
  my $file = shift;

  my $self = {};
  $self->{OUT_FH} = undef;       # Output file handle
  $self->{OUTFILE} = undef;      # Output filename
  $self->{FILE} = $file;         # Current filename processing
  $self->{WEB} = undef;          # Web_File object
  bless( $self );
  return( $self );
}

#
# Methods
#

sub File {
#
# return the current file data
#
  my $self = shift;

  my $file = $self->{FILE};
  return( $file );
}

sub Outfile {
#
# return the outfile data
#
  my $self = shift;

  my $file = $self->{OUTFILE};
  return( $file );
}

sub Next {
#
# Find the next file in the series
#
  my $self = shift;

  my $next_ref = $self->{WEB}->next;
  return( $next_ref );
}

sub Web {
#
# Find web-file object that's within the HTML2Word object
#
  my $self = shift;

  my $web = $self->{WEB};
  return( $web );
}

sub open_out {
#
#  Open the output file and new up a Web_File object to use
#
  my $self = shift;
  my $outname = shift;
  my $outnum = shift;

  my $nm = ref($self) . "::open_out";
  if ( ! defined($outname) ) {
    die "Error:($nm) outname not given to open_out\n";
  }
  if ( ! defined($outnum) ) {
    die "Error:($nm) outnum not given to open_out\n";
  }
  my $outfile = $outname.$outnum.".html";
  $self->{OUTFILE} = $outfile;
  print "Open new output file: $outfile\n";
  open( OUT, ">$outfile" ) || die "Could not open $outfile";
  print OUT "<HTML>\n";
  print OUT '<BODY  LINK="#000000" VLINK="#800080" BGCOLOR="#ffffff">';
  print OUT "\n";
  $self->{OUT_FH} = \*OUT;
  #
  # Also new up a Web_file object for use by the conversion process
  # 
  $self->{WEB} = Web_File->new( );
}

sub close_out {
#
# Close the output file
#
  my $self = shift;

  my $fh = $self->{OUT_FH};
  print $fh "</BODY>\n";
  print $fh "</HTML>\n";
  close( $fh );
  my $outfile = $self->Outfile;
  print "Realign file \n";
  system( "./perl_mod/realign.pl $outfile" );
}

sub strip_file_headers {
#
# Go through the given input file and strip out the navigation buttons
# and other item's that are meant for the Web-based version only.
#
  my $self = shift;
  my $file = shift;

  my $outfile = $self->{OUTFILE};

  my $fh_out = $self->{OUT_FH};
  open( FILE, "<$file" ) || die "Could not open file $file";
  my $nlines = 0;
  my $ntail = 0;
  $self->Web->restart( $file );
  $self->Web->file( $file );
  #
  # Go through all the lines in the file
  #
  while( <FILE> ) {
    #
    # Ignore the header lines
    #
    if ( $nlines == 0 ) {
      if( $self->Web->test_if_begin_file( $_, $fh_out ) ) {
        $nlines++;
        $self->Web->begin_file( $fh_out );
      }
      else {
        #
        # Get the title (HTML can be upper or lower case)
        #
        if ( /\<[tT][iI][tT][lL][eE]>(.+)\<\/[tT][iI][tT][lL][eE]\>/ ) {
          $self->Web->title( $1 );
        }
      }
    }
    elsif ( $ntail > 0 ) {
      $ntail++;
    }
    else {
      #
      # Ignore the tail lines
      #
      if( &Web_File::test_if_end_file( $_ ) ) {
        if ( /(.+)(<!)/ ) {
          print $fh_out "$1\n";
        }
        $ntail++;
      }
      else {
        # Increment # of lines and print to output file
        $nlines++;
        #
        print $fh_out "$_";
      }
    }
  }
  close( FILE);
  &Web_File::end_file( $file, $fh_out );
}
