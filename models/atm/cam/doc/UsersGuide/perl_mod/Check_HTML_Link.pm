#
#       Check_HTML_Link.pm perl5 object
#
#       This module examines links in a HTML file and checks to ensure
#	 that the referenced files and targets exist.
#
#       Version control:
#
#       $Id: Check_HTML_Link.pm,v 1.5 2000/05/02 21:50:58 erik Exp $
#
#       Methods:
#
#	Constructor (new) : 	Creates a Check_HTML_Link object, opens relevant
#				files etc.  The file input to it is the file
#				that is being read.
#	$clink = Check_HTML_Link->new( $file );
#
#	Methods:
#
#	check :		Checks if the input line read from the file is
#				a HTML reference.  If it is a reference it figures
#				out the referenced file and target.  Then it ensures
#				that both the referenced file and target exists.
#	$clink->check( $_ );
#
package Check_HTML_Link;
use strict;

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Constructor
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub new {
#
# Constructor
#
  my $self = {};
  my $name   = shift;
  my $infile = shift;

  my $file = "links.out";
  my $bad_file = "BAD_links.out";
  open( OUT, ">$file" ) || die "Could not open for output: $file";
  open( BAD_OUT, ">$bad_file" ) || die "Could not open for output: $file";

  $self->{FILE} = $infile;       # Input filename examining
  $self->{OUT_FILE} = $file;     # Output filename of list of links
  $self->{OUT} = \*OUT;          # Output filehandle
  $self->{BAD_FILE} = $bad_file; # Output filename of list of bad links
  $self->{BAD_OUT} = \*BAD_OUT;  # Output filehandle to bad links list
  bless( $self );
  return( $self );
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sub file {
#
# Get the list of links filename
#
  my $self = shift;

  my $file = $self->{FILE};
  return( $file );
}

sub out_file {
#
# Get the list of links output filename
#
  my $self = shift;

  my $file = $self->{OUT_FILE};
  return( $file );
}

sub out_fh {
#
# Get the list of bad-links output file-handle
#
  my $self = shift;

  my $out_fh = $self->{OUT};
  return( $out_fh );
}

sub bad_file {
#
# Get the list of bad-links output filename
#
  my $self = shift;

  my $file = $self->{BAD_FILE};
  return( $file );
}

sub bad_out_fh {
#
# Get the list of bad-links output file-handle
#
  my $self = shift;

  my $out_fh = $self->{BAD_OUT};
  return( $out_fh );
}

sub check {
#
# Check a link.  Input a line read in from a HTML file, check to
# see if the line is a HTML link, parse the link out, check that
# the referenced file and target exist.  If not report an error.
# Also write the list of links to a file.
#
  my $self = shift;
  my $line = shift;

  #
  # Check if line has a HTML link reference in it
  #
  my $parsed = undef;
  if ( $line =~ /href[ 	]*=/i ) {
    #
    # Loop over all links in the line
    #
    while( $line =~ /<[Aa][ 	]+[hH][rR][eE][fF][ 	]*=[ 	]*\"*(.+?)\"*[ 	]*>(.*)/i ) {
      $parsed = "True";
      my $file = $1;
      my $target = undef;
      my $link = $2;
      if ( $file =~ /^http\:\/\// ) {
        print "HTTP LINK: Can not check: $file\n";
        print "HTTP LINK: Link: $link\n";
        $file = undef;
        $target = undef;
      }
      elsif ( $file =~ /(.+\.html)\#(.+)/i ) {
        $file = $1;
        $target = $2;
      }
      elsif( $file =~ /^\#(.+)/i ) {
        $target = $1;
        $file = $self->{FILE};
      }
      #
      # if not a HTTP web link go through and check that target exists in file
      #
      if ( defined($file) ) {
        #
        # If referenced file does not exist
        #
        my $bad_out = $self->{BAD_OUT};
        if ( ! -e $file ) {
          print "BAD LINK: $file does not exist!\n";
          print $bad_out "BAD LINK: $file does not exist!\n";
        }
        #
        # Otherwise check if the specified target in the file exists
        #
        else {
          my $check = undef;
          if ( ! defined($target) ) {
            $check = '<[tT][iI][tT][lL][eE]';
          }
          else {
            $check = "[nN][aA][mM][eE][ 	]*=\"*[	 ]*$target";
          }
          my $result = `egrep \'$check\' $file`;
          if ( $result !~ /.+/ ) {
            print "BAD LINK: target $target does not exist in file: $file!\n";
            print $bad_out "BAD LINK: target $target does not exist in file: $file!\n";
          }
          else {
            my $out = $self->{OUT};
            print $out "Link: $file :Target: $target :Link: $link\n";
            print $out "Link: Result: $result\n";
          }
        }
      }
      #
      # Substitute out this link, go to the next link in the line
      #
      $line =~ s/<[aA][ 	]+[hH][rR][eE][fF][     ]*=[   ]*\"*(.+?)\"*[        ]*>//i;
    }
    #
    # If couldn't parse the link out so report an error
    #
    if ( !defined( $parsed ) ) {
      print "Reference line: $line\n";
      print "Possibly the referenced line does not include \<A ?\n";
      die "Found reference to href but could not expand it out!";
    }
  }
  #
  # Return TRUE if a link was found
  #
  if ( defined( $parsed ) ) {
    return( 1 );
  }
  else {
    return( 0 );
  }
}


1 # to keep use or require happy
