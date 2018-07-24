#
#	web_file perl5 object
#
#	This module helps deal with CCM web_files being converted back
#	and forth from MS-Word.  Basically it helps keep track of the
#	internal filenames, and begin and end of files.
#
#	Version control:
#
#	$Id: Web_File.pm,v 1.10.2.1 2002/05/30 17:40:35 erik Exp $
#
#	Methods:
#
#	new ----- Get a new object, assign all values to NULL
#		usage: $web = Web_File->new;
#	find_begin_of_file -- Test if current line is the beginning of
#				an internal file.  If so fill up the object
#				with data on current, next file etc.
#		usage: if ( $web->find_begin_of_file( $line ) )...
#	file ---- Get the current filename.
#		usage: $file = $web->file;
#	next ---- Get the next filename.
#		usage: $file = $web->next;
#	prev ---- Get the previous filename.
#		usage: $file = $web->prev;
#	up ------ Get the top-of-section (up) filename.
#		usage: $file = $web->up;
#	title --- Get the title
#		usage: $title = $web->title;
#	link ---- Get or set the sub-links
#		usage: $link = $web->link; or $web->link( $link );
#      test_if_begin_file --- Test if got beyond the header stuff in file.
#             usage: if ( $web->test_if_begin_file( $_, \*OUT ) );
#	test_end_of_file -- Test if this line is the end of current file.
#             usage: if ( $web->test_end_of_file( $_, \*FILE, \*OUT ) ) ...
#	end_this_file --- Write the ending marks to the file.
#		usage: $web->end_this_file( \*FILE );
#	begin_new_file - Write the beginning marks to the file.
#		usage: $web->begin_new_file( \*FILE );
#	buttons ---- Write the buttons out to the file.
#		usage: $web->buttons( \*FILE );
#	end_file --- Write the end internal file markers to the file
#		usage: &end_file;
#	begin_file --- Write the begin internal file markers to the file
#		usage: &begin_file( );
#	up_file --- Write the up internal file markers to the file
#		usage: &up_file;
#      test_if_end_file --- Test if got to end of body of the file
#             usage: if ( &Web_File::test_if_end_file( $_ ) );
#	find_headers -- Find headers and markers to use in the table of contents.
#		(for pages that contain more than one header).
#		usage: $web->find_headers( $_ );
#
use strict;
require "Section.pm";
package Web_File;

#-----------------------------------------------------------
# Constructor
#-----------------------------------------------------------
sub new 
#
# Construct a web_file object
#
# usage: $web = Web_File->new;
#
{
  my $self = {};

  print "new up a Web_File object: \n";
  $self->{FILE} = undef;  # filename of current internal file
  $self->{NEXT} = undef;  # filename of next internal file
  $self->{PREV} = undef;  # filename of previous internal file
  $self->{TITLE}= undef;  # title of current internal file
  $self->{UP}   = undef;  # filename of head of this section
  $self->{LINK} = undef;  # List of sub-header links or "TOC" if table of contents
  $self->{MARKER} = undef;# Last marker (target) read in the page
  $self->{PRE_UP} = undef;# filename of the up section for previous file
  $self->{PRE_TITLE} = undef;    # title of the previous file
  $self->{SECTION} = Section->new( );# names and such for where we are at in the doc.
  $self->{TOC} = Section->new( );# table of contents links.
  $self->{TOC_SET} = undef;      # if table of contents has been set.
  $self->{START}= 0;             # if got to read the starting file or not
  bless( $self );
  return( $self );
}

#-----------------------------------------------------------
# Methods
#-----------------------------------------------------------

sub display {
#
# Display the contents of the object
#
  my $self = shift;

  my $file = $self->{FILE};
  my $next = $self->{NEXT};
  my $prev = $self->{PREV};
  my $title = $self->{TITLE};
  my $up = $self->{UP};
  my $link = $self->{LINK};
  my $pre_up = $self->{PRE_UP};
  my $pre_title = $self->{PRE_TITLE};
  my $TOC_SET = $self->{TOC_SET};
  my $section = $self->{SECTION};
  my $TOC     = $self->{TOC};
  my $start = $self->{START};
  my $marker = $self->{MARKER};

  if( defined($file) ) {
    print "File      = $file\n";
  }
  if( defined($next) ) {
    print "Next      = $next\n";
  }
  if( defined($prev) ) {
    print "Prev      = $prev\n";
  }
  if( defined($up) ) {
    print "Up        = $up\n";
  }
  if( defined($link) ) {
    print "Link      = $link\n";
  }
  if( defined($pre_up) ) {
    print "pre_up    = $pre_up\n";
  }
  if( defined($marker) ) {
    print "marker    = $marker\n";
  }
  if( defined($pre_title) ) {
    print "pre_title = $pre_title\n";
  }
  if( defined($TOC_SET) ) {
    print "TOC_SET   = $TOC_SET\n";
  }
  if( defined($start) ) {
    print "start     = $start\n";
  }
  if( defined($section) ) {
    print "Section: \n";
    $section->display;
  }
  if( defined($TOC) ) {
    print "TOC: \n";
    $TOC->display;
  }
}

sub restart
#
# Reset the start message to zero when we go to a new file
#
{
  my $self = shift;
  my $file = shift;

  #
  # Set as if we had just ended this
  #
  $self->{PREV} = $file;
  $self->{PRE_TITLE}= undef;
  $self->{PRE_UP}   = undef;
  #
  # Reset internal data to undefined
  #
  $self->{FILE} = undef;  # filename of current internal file
  $self->{NEXT} = undef;  # filename of next internal file
  $self->{TITLE}= undef;  # title of current internal file
  $self->{UP}   = undef;  # filename of head of this section
  $self->{START}= 0;      # if got to read the starting file or not
  return;
}
#------------------------------------------------------------------
# Extracting data from the object
#------------------------------------------------------------------

sub file
#
# Get the filename of the current file
#
# usage: $data = $web->file;
#
{
  my $self = shift;
  my $file = shift;

  if ( defined( $file ) ) {
    $self->{FILE} = $file;
    return;
  }
  else {
    $file = $self->{FILE};
    return( $file );
  }
}

sub marker
#
# Get the marker of the current file
#
{
  my $self = shift;

  my $marker = $self->{MARKER};
  return( $marker );
}

sub prev
#
# Get the filename of the previous file
#
# usage: $data = $web->prev;
#
{
  my $self = shift;

  my $file = $self->{PREV};
  return( $file );
}

sub next
#
# Get the filename of the next file
#
# usage: $data = $web->next;
#
{
  my $self = shift;
  my $next = shift;

  if ( defined( $next ) ) {
    $self->{NEXT} = $next;
    return;
  }
  else {
    $next = $self->{NEXT};
    return( $next );
  }
}

sub Table_Of_Cont
#
# Get the TOC (table of contents section) of the current file
#
# usage: $data = $web->Table_Of_Cont;
#
{
  my $self = shift;

  my $file = $self->{TOC};
  return( $file );
}

sub Section
#
# Get the Section of the current file
#
# usage: $data = $web->Section;
#
{
  my $self = shift;

  my $file = $self->{SECTION};
  return( $file );
}

sub Set_TOC
#
# Set the TOC flag so additional links won't be added after the
# TOC has been calculated.
#
{
  my $self = shift;

  $self->{TOC_SET} = "True";
}

sub up
#
# Get the filename of the up (top of section) file
#
# usage: $data = $web->up;
#
{
  my $self = shift;

  my $file = $self->{UP};
  return( $file );
}

sub title
#
# Get the title of the page
#
# usage: $data = $web->title;
#
{
  my $self = shift;
  my $title= shift;

  if ( defined( $title) ) {
    $self->{TITLE} = $title;
    return;
  }
  else {
    $title= $self->{TITLE};
    return( $title);
  }
}


sub link
#
#	Get or set the sub-links
# usage: $link = $web->link; or $web->link( $link );
#
{
  my $self = shift;
  my $link= shift;

  if ( defined( $link) ) {
    $self->{LINK} = $link;
    return;
  }
  else {
    $link= $self->{LINK};
    return( $link);
  }
}

#------------------------------------------------------------------
# General methods...
#------------------------------------------------------------------

sub test_if_begin_file
#
# When reading separate files check if got beyond the header stuff yet
# Also get the up-file reference
#
{
  my $self = shift;
  my $line = shift;
  my $fh = shift;
#Sample:<A HREF="UG-8.html"><IMG ALIGN = BOTTOM SRC = "images/up.gif"></A>
#  Get the up-file reference
#
  if ( $line =~ /up[g]*\.gif/ ) {
#   Do a non-greedy case-insensitive match that has the nearest up.gif in it
#   HTML can be upper or lower case
    if ( /^.*\<[aA]\s+[hH][rR][eE][fF]\s*="(.+?)"\>(.+?)up[g]*\.gif/i ) {
      $self->{UP} = $1;
    }
    else {
      print "UP-file reference: $_";
      die "did not get the up reference:";
    }
  }
#Sample:<A HREF="UG-8.html"><IMG ALIGN = BOTTOM SRC = "images/next.gif"></A>
#  Get the next-file reference
#
  if ( $line =~ /next[g]*\.gif/ ) {
#   Do a non-greedy case-insensitive match that has the nearest next.gif in it
#   HTML can be upper or lower case
    if ( /^.*\<[aA]\s+[hH][rR][eE][fF]\s*="(.+?)"\>(.+?)next[g]*\.gif/i ) {
      $self->{NEXT} = $1;
    }
    else {
      print "NEXT-file reference: $_";
      die "did not get the next reference:";
    }
  }
  if ( $line =~ /<!Beginning_of_the_page:/ ) {
    return( 1 );  # found beginning of the page
  }
  else {
    return( 0 );  # haven't found it yet...
  }
}

sub make_sub_links
#
# Mark the sub-section links at the end of the file
#
# usage: $web->make_sub_links( \*OUT )
#
{
  my $self = shift;
  my $out_fh = shift;

  print $out_fh "<HR>\n";
  print $out_fh "<CENTER><H2>Sub Sections</H2></CENTER>\n";
  print $out_fh "<DL>\n";
  my $link = $self->{LINK};
  while( $link->File !~ /NULL/ ) {
    $link->write_HTML_node( $out_fh );
    $link = $link->Sub_Section;
  }
  print $out_fh "</DL>\n";
  print $out_fh "<P>\n";
  $self->{LINK}->DESTROY;
  $self->{LINK} = undef;
  return;
}

sub test_end_of_file
#
# Test for the end of the current file
#
# usage: if ( $web->test_end_of_file( $_, \*FILE, \*OUT ) ) ...
#
{
  my $self = shift;
  my $line = shift;
  my $in_fh = shift;
  my $out_fh = shift;

  if ( $line =~ /\<\!end_of_file:/ ) {
    my $file = $self->file;
    #
    # Close out file
    #
    if ( $file =~ /table_of_contents/ ) {
      $self->{LINK} = "TOC";
    }
    my $next = $self->next;
    my $set = undef;
    #
    # If reached a file that has no "next" set the Table of Contents
    #
    if ( ! defined( $next ) ) {
      $set = 1;
    }
    $self->end_this_file( $line, $in_fh, $out_fh );
    #
    # Actually set the table of contents now that current file
    # has been added.
    #
    if ( defined($set) ) {
      print "Set the Table of Contents\n";
      $self->Set_TOC;
    }
    return( 1 );
  }
  else {
    return( 0 );
  }
  return;
}


sub find_begin_of_file
#
# Find the beginning of this file
#
# usage: if ( $web->find_begin_of_file( $_ ) ) ...
#
{
  my $self = shift;
  my $line = shift;

  my $file = undef;

  if ( $line =~ /\<\!begin_of_file: (.+)>$/ ) {
    $self->{FILE} = $1;
    $self->{START} = 1;
  }
  elsif ( $self->{START} == 1 ) {
    my $set = 0;
    if ( /\<\!next_file: (.+)>$/ ) {
      $self->{NEXT} = $1;
      $set = 1;
    }
    if ( $line =~ /\<\!up_file: (.+)>$/ ) {
      $self->{UP} = $1;
      $set = 1;
    }
    if ( $line =~ /\<\!title_of_file: (.+)>$/ ) {
      $self->{TITLE} = $1;
      $set = 1;
    }
    if ( $set != 1 ) {
      print "File: $self->{FILE}\n";
      return( 1 );
    }
  }
  return( 0 );
}

sub buttons
#
# Draw the navigation buttons on the page
#
# usage: $web->buttons( \*FILE );
#
{
  my $self = shift;
  my $fh   = shift;

  my $file = $self->file;
  my $next = $self->next;
  my $prev = $self->prev;
  my $up   = $self->up;
  my $cont  = 0;
  my $pr    = 0;
  my $nx    = 0;
  my $top   = 0;
  my $secup = 0;
  if ( $file =~ /index/ ) {
    $cont  = 1;
    $nx    = 1;
  }
  elsif ( $file =~ /index/ ) {
    $self->{PREV} = "search.html";
    $prev  = $self->{PREV};
    $pr    = 1;
    $cont  = 1;
    $nx    = 1;
  }
  elsif ( $file =~ /table_of_contents/ ) {
    $pr    = 1;
    $nx    = 1;
    $top   = 1;
  }
  elsif ( defined( $up ) ) {
    $cont  = 1;
    $pr    = 1;
    $nx    = 1;
    $top   = 1;
    $secup = 1;
  }
  else {
    $cont  = 1;
    $pr    = 1;
    $nx    = 1;
    $top   = 1;
  }
  if ( ! defined( $next ) ) {
    $nx    = 0;
  }
  print $fh "<BR>\n";
  #
  # Add the buttons in one by one if they are referenced above...
  #
  if ( $nx == 1 ) {
    print $fh qq(<A HREF="$next"><IMG SRC="images/next.gif" \n);
    print $fh qq(ALT = "Go to next page"\n);
    print $fh qq( ALIGN=BOTTOM></A>\n);
  }
  if ( $pr == 1 ) {
    print $fh qq(<A HREF="$prev"><IMG SRC="images/prev.gif" \n);
    print $fh qq(ALT = "Go to previous page"\n);
    print $fh qq( ALIGN=BOTTOM></A>\n);
  }
  if ( $secup == 1 ) {
    print $fh qq(<A HREF="$up"><IMG SRC="images/up.gif" \n);
    print $fh qq(ALT = "Go to top of this section"\n);
    print $fh qq( ALIGN=BOTTOM></A>\n);
  }
  if ( $top == 1 ) {
    print $fh qq(<A HREF="index.shtml"><IMG SRC="images/top.gif" \n);
    print $fh qq(ALT = "Go to top page"\n);
    print $fh qq( ALIGN=BOTTOM></A>\n);
  }
  if ( $cont == 1 ) {
    print $fh qq(<A HREF="table_of_contents.html"><IMG SRC="images/content.gif" \n);
    print $fh qq(ALT = "Go to table of contents"\n);
    print $fh qq( ALIGN=BOTTOM></A>\n);
  }
  print $fh "<BR>\n";
  return;
}
sub make_section_links
#
# Make the links to the previous sections at the top of the page
#
{
  my $self = shift;
  my $fh   = shift;

  my $section = $self->{SECTION};
  my $toc = $self->{TOC};
  my $file = $self->{FILE};
  my $title= $self->{TITLE};
  my $up   = $self->up;
  if ( $file !~ /index\.html/ ) {
    $section->Make_Section_Link( $fh, $title, $file, $up );
    if ( ! defined($self->{TOC_SET}) ) {
      $toc->Make_TOC_Link( $title, $file, $up );
    }
  }
}

sub begin_new_file
#
# Write out standard opening information to beginning of file
#
# usage: $web->begin_file( \*FILE );
#
{
  my $self = shift;
  my $fh   = shift;

  my $file = $self->file;
  my $title = $self->title;
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Routines to write out in the ordinal file the internal file markers...
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub begin_file
#
#  Mark the beginning of the internal file
#  with a HTML comment
#
{
  my $self = shift;
  my $fh = shift;

  my $file  = $self->{FILE};
  my $next  = $self->{NEXT};
  my $title = $self->{TITLE};
  my $up    = $self->{UP};
  print $fh  "<\!begin_of_file: $file>\n";
  if ( defined($title) ) {
    print $fh  "<\!title_of_file: $title>\n";
  } else {
    die "Title not defined";
  }
  if ( defined($up) ) {
    my $up_file = $up;
    print $fh "<\!up_file: $up_file>\n";
  }
  if ( defined($next) ) {
    my $next_file = $next;
    print $fh "<\!next_file: $next_file>\n";
  }
  return;
}

sub up_file
#
#  Mark the file that marks the beginning of the section
#  with a HTML comment
#
{
  my $up_file = shift;
  my $fh = shift;

  print $fh  "<\!up_file: $up_file>\n";
  return;
}

sub end_file
#
#  Mark the end of the internal file
#  with a HTML comment
#
{
  my $file = shift;
  my $fh = shift;

  print $fh  "<\!end_of_file: $file>\n";
  return;
}

sub test_if_end_file
#
# When reading separate files check if got beyond the body
#
{
  my $line = shift;

  if ( $line =~ /<!End_of_the_page:/ ) {
    return( 1 );  # found end of the page
  }
  else {
    return( 0 );  # haven't found it yet...
  }
}

sub get_links
#
# Get the list of links from the end of the given file
#
{
  my $in_fh = shift;
  my $out_fh = shift;

}

sub find_headers
#
# Find headers and markers to use in the table of contents.
# (for pages that contain more than one header).
#
{
  my $self = shift;
  my $line = shift;

  # HTML can be upper or lower case
  if ( $line =~ /<a\s+name\s*=\s*"*([^>]+?)"*>/i ) {
    $self->{MARKER} = $1;
  }
  my $new_line = $line;
  $new_line =~ s/<a\s+name\s*=\s*"*([^>]+?)"*>([^<]*?)<\/a>//gi; 
  my $header;
  if ( defined($2) ) {
    $header = $2;
  } else {
    $header = "";
  }
  if ( $new_line =~ /<h[1-9]>(.*?)<\/h/i ) {
    if ( defined($1) ) {
      $header = "$header$1";
    }
    my $title = $self->{TITLE};
    if ( (! defined($self->{TOC_SET})) && ($header !~ /$title/) ) {
      my $marker = $self->{MARKER};
      my $file = $self->file;
      my $target = undef;
      if ( defined( $marker ) ) {
        $target = "$file#$marker";
      }
      else {
        $target = "$file";
      }
      my $toc = $self->Table_Of_Cont;
      my $up = $marker;
      $toc->Make_TOC_Link( $header, $target, $up );
    }
  } elsif ( $new_line =~ /<h[1-9]>/i ) {
     my $file = $self->file;
     print "Can not parse line: $new_line, file: $file";
     die "Trouble parsing line, headers need to be on one line";
  }
  return;
}



1   # to make use or require happy
