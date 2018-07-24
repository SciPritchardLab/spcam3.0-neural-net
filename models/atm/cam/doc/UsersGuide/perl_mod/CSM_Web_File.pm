#
#	CSM_Web_File.pm perl5 object
#
#	This module helps deal with CSM web_files being converted back
#	and forth from MS-Word.  Basically it helps keep track of the
#	internal filenames, and begin and end of files.
#
#	Version control:
#
#	$Id: CSM_Web_File.pm,v 1.3 1999/03/25 22:36:35 erik Exp $
#
#	Methods:
#
#	Inherient all methods from the Web_file object other than
#	begin_new_file that writes the header information for a file, and
#	end_this_file that writes the ending information.
#
#
package CSM_Web_File;
@ISA = "Web_File";
use Web_File;

#-----------------------------------------------------------
# Constructor
#-----------------------------------------------------------
#
# Constructor
#
  sub new {
    my $class = shift;
    my $self = {};
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
  print $fh "<!---  csm_doc/$file					----->";
  print $fh "\n";
  print $fh '<!---								----->';
  print $fh "\n";
  print $fh '<!---  HTML version of the documentation on the NCAR CSM	----->';
  print $fh "\n";
  print $fh '<!---  (Climate System Model) version 2.0	----->';
  print $fh "\n";
  print $fh '<!---								----->';
  print $fh "\n";
  print $fh '<!---  Version control information:				----->';
  print $fh "\n";
  print $fh '<!---								----->';
  print $fh "\n";
  print $fh '<!---  $Id: CSM_Web_File.pm,v 1.3 1999/03/25 22:36:35 erik Exp $			----->';
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
  print $fh qq(<META NAME="description" CONTENT="$title">\n);
  print $fh qq(<META NAME="resource-type" CONTENT="document">\n);
  print $fh qq(<META NAME="distribution" CONTENT="global">\n);
  print $fh qq(<META HTTP-EQUIV="Content-Type" CONTENT="text/html; charset=iso_8859_1">\n);
  print $fh qq(<LINK REL="STYLESHEET" HREF="csm_users_guide.css">\n);
  my $next = $self->next;
  my $prev = $self->prev;
  my $up   = $self->up;
  if ( defined( $up ) ) {
    print $fh qq(<LINK REL="up" HREF="$up">\n);
  }
  if ( defined( $prev ) ) {
    print $fh qq(<LINK REL="previous" HREF="$prev">\n);
  }
  if ( defined( $next ) ) {
    print $fh qq(<LINK REL="next" HREF="$next">\n);
  }
  print $fh '</HEAD>';
  print $fh "\n";
  print $fh qq(<BODY  TEXT="#990000" BACKGROUND="images/Background.gif" LINK="#000099" );
  print $fh "\n";
  print $fh qq(  VLINK="#000099" ALINK="#FF0000">);
  print $fh "\n";
  print $fh "<A NAME=TOP_OF_PAGE><P>\n";
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
  #
  # Add comment that this is the beginning of the page
  #
  print $fh '<!the navigation buttons and general format are edited>';
  print $fh "\n";
  print $fh '<!in the word2html.pl and Web_File.pm Perl5 script.>';
  print $fh "\n";
  print $fh '<!Beginning_of_the_page: -- do not edit anything above!!!>';
  print $fh "\n";
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
  print $fh qq(&nbsp;<A HREF="index.html"><IMG SRC="images/CSMLogo.gif" \n); 
  print $fh qq(ALT = "Search for keywords in the CSM2.0 Users Guide" \n);
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
  print $fh '<A href="mailto:csm@ucar.edu">csm@ucar.edu</A> .';
  print $fh "\n";
  print $fh '<HR>';
  print $fh "\n";
  print $fh "<!----------------------------------------------------------------------------->\n";
  print $fh "<!---                        last update message                            --->\n";
  print $fh "<!----------------------------------------------------------------------------->\n";
  print $fh "<BR><BR> <HR>\n";
  print $fh "<CENTER>\n";
  print $fh qq(<script language="JavaScript">\n);
  print $fh " <!--hide script from old browsers\n";
  print $fh qq(  document.write("Updated: " + document.lastModified );\n);
  print $fh "   // end hiding -->\n";
  print $fh "</script>\n";
  print $fh "</CENTER>\n";
  print $fh "<HR><BR><BR>\n";

  print $fh '<CENTER><ADDRESS>$Name:  $ $Revision: 1.3 $ $Author: erik $</ADDRESS></CENTER>';
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
  #
  # Set the previous file to the current file
  #
  $self->{PREV} = $self->{FILE};
  $self->{PRE_TITLE}= $self->{TITLE};  # title of previous internal file
  $self->{PRE_UP}   = $self->{UP}; # up-filename of previous internal file
  #
  # Reset internal data to undefined
  #
  $self->{FILE} = undef;  # filename of current internal file
  $self->{NEXT} = undef;  # filename of next internal file
  $self->{TITLE}= undef;  # title of current internal file
  $self->{UP}   = undef;  # filename of head of this section
  $self->{START}= 0;  # if got to read the starting file or not
  return;
}



1   # to make use or require happy
