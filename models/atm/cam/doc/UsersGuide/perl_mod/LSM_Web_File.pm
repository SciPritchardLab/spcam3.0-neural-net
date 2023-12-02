#
#	LSM_Web_File.pm perl5 object
#
#	This module helps deal with LSM web_files being converted back
#	and forth from MS-Word.  Basically it helps keep track of the
#	internal filenames, and begin and end of files.
#
#	Version control:
#
#	$Id: LSM_Web_File.pm,v 1.11 2001/11/01 16:14:55 erik Exp $
#
#	Methods:
#
#	Inherient all methods from the Web_file object other than
#	begin_new_file that writes the header information for a file, and
#	end_this_file that writes the ending information.
#
#
package LSM_Web_File;
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
  print $fh qq(<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2//EN">\n);
  print $fh '<!------------------------------------------------------------------->';
  print $fh "\n";
  print $fh '<!---								----->';
  print $fh "\n";
  print $fh "<!---  lsm_doc/$file					----->";
  print $fh "\n";
  print $fh '<!---								----->';
  print $fh "\n";
  print $fh '<!---  HTML version of the documentation on the NCAR LSM	----->';
  print $fh "\n";
  print $fh '<!---  land surface model version 1.1 as coupled to the CCM	----->';
  print $fh '<!---  atmospheric model.					----->';
  print $fh "\n";
  print $fh '<!---								----->';
  print $fh "\n";
  print $fh '<!---  Version control information:				----->';
  print $fh "\n";
  print $fh '<!---								----->';
  print $fh "\n";
  print $fh '<!---  $Id: LSM_Web_File.pm,v 1.11 2001/11/01 16:14:55 erik Exp $			----->';
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
  print $fh qq(<META NAME="keywords" CONTENT="lsm_users_guide">\n);
  print $fh qq(<META NAME="resource-type" CONTENT="document">\n);
  print $fh qq(<META NAME="distribution" CONTENT="global">\n);
  print $fh qq(<META HTTP-EQUIV="Content-Type" CONTENT="text/html; charset=iso_8859_1">\n);
  print $fh qq(<LINK REL="STYLESHEET" HREF="lsm_users_guide.css">\n);
  my $next = $self->next;
  my $prev = $self->prev;
  my $up   = $self->up;
  print $fh qq(<LINK REL="previous" HREF="table_of_contents.html">\n);
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
  print $fh '<BODY BGCOLOR = "WHITE">';
  print $fh "\n";
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
  print $fh '<P>';
  print $fh "\n";
  print $fh '<HR>';
  print $fh "\n";
  $self->buttons( $fh );
  #
  # Put the globe and search button at the bottom of the page...
  #
  print $fh '<P>';
  print $fh "\n";
  print $fh qq(&nbsp;<A HREF="index.shtml"><IMG SRC="images/globe.gif" \n); 
  print $fh qq(ALT = "Search for keywords in the LSM1.1 as coupled to CAM1.8 Users Guide" \n);
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
  print $fh '<A href="mailto:mvertens@ucar.edu">mvertens@ucar.edu</A> .';
  print $fh "\n";
  print $fh '<HR>';
  print $fh "\n";
  print $fh '<ADDRESS>$Name:  $ $Revision: 1.11 $ $Date: 2001/11/01 16:14:55 $ $Author: erik $</ADDRESS>';
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
