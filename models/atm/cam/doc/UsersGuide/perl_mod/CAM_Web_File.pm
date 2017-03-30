#
#	CAM_Web_File.pm perl5 object
#
#	This module helps deal with CAM web_files being converted back
#	and forth from MS-Word.  Basically it helps keep track of the
#	internal filenames, and begin and end of files.
#
#	Version control:
#
#	$Id: CAM_Web_File.pm,v 1.1.2.2 2002/05/30 17:40:33 erik Exp $
#
#	Methods:
#
#	Inherient all methods from the Web_file object other than
#	begin_new_file that writes the header information for a file, and
#	end_this_file that writes the ending information.
#
#
use strict;
package CAM_Web_File;
@CAM_Web_File::ISA = "Web_File";
use Web_File;

#-----------------------------------------------------------
# Constructor
#-----------------------------------------------------------
#
# Constructor
#
  sub new {
    my $class = shift;
    my $tag = shift;
 
    my $self = $class->SUPER::new( );
    bless( $self, $class );
    return( $self );
  }

#-----------------------------------------------------------
# Methods
#-----------------------------------------------------------

sub begin_new_file
#
# Write out standard opening information to beginning of file
#
# usage: $web->begin_new_file( \*FILE );
#
{
  my $self = shift;
  my $fh   = shift;

  my $file = $self->file;
  my $title = $self->title;
  print $fh '<!------------------------------------------------------------------->';
  print $fh "\n";
  print $fh '<!---								----->';
  print $fh "\n";
  print $fh "<!---  cam_doc/$file					----->";
  print $fh "\n";
  print $fh '<!---								----->';
  print $fh "\n";
  print $fh '<!---  HTML version of the documentation on the NCAR global	----->';
  print $fh "\n";
  print $fh '<!---  atmospheric model CAM2.0.					----->';
  print $fh "\n";
  print $fh '<!---								----->';
  print $fh "\n";
  print $fh '<!---  Version control information:				----->';
  print $fh "\n";
  print $fh '<!---								----->';
  print $fh "\n";
  print $fh '<!---  $Id: CAM_Web_File.pm,v 1.1.2.2 2002/05/30 17:40:33 erik Exp $			----->';
  print $fh "\n";
  print $fh '<!---								----->';
  print $fh "\n";
  print $fh '<!------------------------------------------------------------------->';
  print $fh "\n";

  print $fh '<HTML>';
  print $fh "\n";
  print $fh '<HEAD>';
  print $fh "\n";
  print $fh "<TITLE>$title</TITLE>";
  print $fh "\n";
  print $fh '</HEAD>';
  print $fh "\n";
  print $fh '<BODY BGCOLOR = "WHITE">';
  print $fh "\n";
  print $fh '<A NAME=TOP_OF_PAGE>';
  print $fh '<P>';
  print $fh "\n";
  print $fh '<HR>';
  print $fh "\n";
  print $fh qq(&nbsp;<A HREF="#BOTTOM_OF_PAGE"><IMG SRC="images/bottom_of_page.gif" \n); 
  print $fh qq(ALT = "Go to the bottom of this page. See the search engine and sub-section links."\n);
  print $fh qq(ALIGN=BOTTOM></A>\n);
  $self->buttons( $fh );
  $self->make_section_links( $fh );
  print $fh '<HR>';
  print $fh "\n";
  $self->SUPER::begin_new_file( $fh );
  return;
}

sub end_this_file
#
#  Write out closing information to end of file
#
# usage: $web->end_this_file( $_, \*INFILE, \*FILE );
#
{
  my $self = shift;
  my $line    = shift;
  my $in_fh   = shift;
  my $fh   = shift;

  my $file = $self->{FILE};
  #
  # Add comment that this is the end of the page
  #
  print $fh '<!End_of_the_page: -- do not edit anything below!!!>';
  print $fh "\n";
  if ( defined($self->{LINK}) ) {
    if ( $self->{LINK} =~ /TOC/ ) {
      #
      # Write out table of contents links
      #  
      my $toc = $self->Table_Of_Cont;
      my $toc_flag = "True";
      $toc->write_HTML( $fh, $toc_flag  );
      $self->{LINK} = undef;
    }
    else {
      #
      # Make sub-section links
      #  
      $self->make_sub_links( $fh );
    }
  }
  #
  # Add the navigation buttons to the bottom of the page...
  #
  print $fh '<A NAME=BOTTOM_OF_PAGE>';
  print $fh '<P>';
  print $fh "\n";
  print $fh '<HR>';
  print $fh "\n";
  print $fh qq(&nbsp;<A HREF="#TOP_OF_PAGE"><IMG SRC="images/top_of_page.gif" \n); 
  print $fh qq(ALT = "Go to the top of this page. See links to previous section headers."\n);
  print $fh qq(ALIGN=BOTTOM></A>\n);

  $self->buttons( $fh );
  #
  # Put the globe and search button at the bottom of the page...
  #
  print $fh '<P>';
  print $fh "\n";
  print $fh qq(&nbsp;<A HREF="search.html"><IMG SRC="images/cam.jpg" \n); 
  print $fh qq(ALT = "Search for keywords in the CAM2.0 Users Guide"\n);
  print $fh qq(ALIGN=BOTTOM>Search page</A>\n);
  print $fh "\n";
  print $fh '<P>';
  #
  # Add the mail interface, and revision/dates
  #
  print $fh "\n";
  print $fh '<P>';
  print $fh "\n";
  print $fh "Questions on these pages can be sent to...\n";
  print $fh '<A href="mailto:erik@ucar.edu">erik@ucar.edu</A> .';
  print $fh "\n";
  print $fh '<HR>';
  print $fh "\n";
  print $fh '<ADDRESS>$Name:  $ $Revision: 1.1.2.2 $ $Date: 2002/05/30 17:40:33 $ $Author: erik $</ADDRESS>';
  print $fh "\n";
  #
  # close out the page with ending HTML flags
  #
  print $fh '<HR\>';
  print $fh "\n";
  print $fh '</BODY>';
  print $fh "\n";
  print $fh '</HTML>';
  print $fh "\n";
  close( $fh );
  $self->SUPER::end_this_file( $line, $in_fh, $fh );
  return;
}



1   # to make use or require happy
