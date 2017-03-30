#
#	Section.pm
#	Erik Kluzek
#	Mar/3/1998
#
#	Perl5 code to deal with a section object tracking filenames and 
#	titles for each section.
#
#	$Id: Section.pm,v 1.9.2.1 2002/05/30 17:40:34 erik Exp $
#

package Section;
use strict;

sub new {
#
# Section constructor object
#
  my $self = {};
  $self->{TITLE} = undef;       # Title of this section.
  $self->{FILE} = "NULL";       # Filename that goes with this section
  $self->{UP}   = undef;        # Filename to current up file reference
  $self->{LEVEL} = undef;       # Header level that this section is at
  $self->{LAST_LEVEL} = undef;  # Level of previous section
  $self->{LAST_UP} = undef;     # up-file reference of previous section
  $self->{SUB_SECTION} = undef; # Next sub-section link to another Section object
  bless( $self );
  return( $self );
}
#------------------------------------------------------------------------------
# Methods
#------------------------------------------------------------------------------

sub File {
#
# Get the filename of this section
#
  my $self = shift;

  my $file = $self->{FILE};
  return( $file );
}

sub File_Name {
#
# Get the filename (without the marker) of this section
#
  my $self = shift;

  my $file = $self->{FILE};
  if ( $file =~ /(.+?)\#(.+?)/ ) {
    $file = $1;
  }
  return( $file );
}

sub Title {
#
# Get the title of this section
#
  my $self = shift;

  my $file = $self->{TITLE};
  return( $file );
}

sub Level {
#
# Get the level (hierarchy) of this section
#
  my $self = shift;

  my $file = $self->{LEVEL};
  return( $file );
}


sub Sub_Section {
#
# Get the next sub-section
#
  my $self = shift;

  my $file = $self->{SUB_SECTION};
  return( $file );
}

sub Title_Level {
#
# Return the level based on the title alone
#
  my $Title = shift;
  my $level = undef;

  if ( ! defined($Title) ) {
    die "Title not entered";
  }
  if ( $Title =~ /^[1234567890]+\.[1234567890]+\.[1234567890]+\.[1234567890]+\.[1234567890]+ / ) {
    $level = 5;
  }
  elsif ( $Title =~ /^[1234567890]+\.[1234567890]+\.[1234567890]+\.[1234567890]+ / ) {
    $level = 4;
  }
  elsif ( $Title =~ /^[1234567890]+\.[1234567890]+\.[1234567890]+ / ) {
    $level = 3;
  }
  elsif ( $Title =~ /^[A-Z1234567890]+\.[1234567890]+ / ) {
    $level = 2;
  }
  elsif ( $Title =~ /^[1234567890]+\. / ) {
    $level = 1;
  }
  else {
    $level = 1;
  }
  return( $level );
}

sub Fill_Section {
#
# Fill the section with the input data
#
  my $self = shift;
  my $Title = shift;
  my $file = shift;
  my $up = shift;

  my $level = undef;
  #
  # If first section, fill the current section
  #
  if ( $self->File =~ /NULL/ ) {
    $level = 0;
    $self->{TITLE} = $Title;
    $self->{FILE} = $file;
    $self->{LEVEL} = 0;
    $self->{SUB_SECTION} = Section->new( );
  }
  #
  # Otherwise fill the sub-section to current section
  #
  else {
    #
    # Get the references to the previous section here
    #
    my $last_up = $self->{UP};
    my $last_level = $self->{LEVEL};
    my $sub = $self->Sub_Section;
    if ( defined($sub->Sub_Section) ) {
      $sub->Sub_Section->DESTROY;
    }
    $sub->{LAST_LEVEL} = $last_level;
    $sub->{LAST_UP} = $last_up;
    $sub->{SUB_SECTION} = Section->new( );
    $sub->{TITLE} = $Title;
    $sub->{FILE} = $file;
    $sub->{UP}   = $up;
    #
    # If a new section LEVEL is set to one
    #
    if ( ! defined( $up ) ) {
      $level = 1;
    }
    else {
      if ( ! defined( $last_up ) ) {
        $level = 2;
      }
      elsif ( $up =~ /$last_up/ ) {
        $level = $last_level;
      }
      else {
        $level = Title_Level( $Title );
      }
    }
    $sub->{LEVEL} = $level;
  }
  return;
}

sub DESTROY {
#
# Mark this object as NULL and delete all it's subsections.
#
  my $self = shift;

  my $sub = $self->{SUB_SECTION};
  if ( defined($sub) ) {
    $sub->DESTROY;
  }
}

sub display_node {
#
# Display the contents of the current node
#
  my $self = shift;

  my $file  = $self->File;
  my $title = $self->Title;
  my $sub   = $self->Sub_Section;
  my $level = $self->Level;
  my $last_lev = $self->{LAST_LEVEL};
  my $up = $self->{UP};
  my $last_up = $self->{LAST_UP};
  print "self:      $self \n";
  print "File:      $file \n";
  print "Title:     $title \n";
  print "Up:        $up \n";
  print "Last Up:   $last_up \n";
  print "Level:     $level \n";
  print "Last Lev:  $last_lev \n";
  print "sub:       $sub \n";
}

sub write_HTML_node {
#
# Write out the HTML version of this node
#
  my $self = shift;
  my $fh_out = shift;

  my $file  = $self->File;
  my $title = $self->Title;
  my $level = $self->Level;
  if ( defined($title) ) {
    if ( $level <= 2 ) {
      print $fh_out "<BR>\n";
    }
    if ( $level == 0 ) {
      print $fh_out "<FONT=+1>\n";
    }
    print $fh_out "<DT>\n";
    my $i = undef;
    for( $i = 1; $i < $level; $i++ ) {
      print $fh_out "&nbsp;&nbsp;&nbsp;\n";
    }
    print $fh_out "<A HREF=$file>$title</A>\n";
    print $fh_out "<DD>\n";
    if ( $level == 0 ) {
      print $fh_out "</FONT>\n";
      print $fh_out "<BR>\n";
    }
  }
}

sub write_HTML {
#
# Write a HTML version of the tree
# Write everything but the last node if not a table of contents.
# If it is a table_of_contents write don't write out section title
# and do write out all nodes.
#
  my $self = shift;
  my $fh_out = shift;
  my $toc = shift;

  my $sub = $self->Sub_Section;
  if( $sub->File =~ /NULL/ ) {
    return;
  }
  #
  # If not table of contents write out the previous section header
  #
  if ( $toc !~ /True/ ) {
    print $fh_out "<CENTER><H3>Previous Section Headers</H3></CENTER>\n";
  }
  print $fh_out "<H4><DL>\n";
  $self->write_HTML_node( $fh_out );
  while( $sub->File !~ /NULL/ ) {
    #
    # If not table of contents write everything but current section
    # If table of contents write all sections
    #
    if ( ($toc =~ /True/) || ($sub->Sub_Section->File !~ /NULL/) ) {
      $sub->write_HTML_node( $fh_out );
    }
    $sub = $sub->Sub_Section;
  }
  print $fh_out "</H4></DL>\n";
}

sub display {
#
# Display the contents of the tree
#
  my $self = shift;

  print "START DISPLAY OF NODES:\n\n";
  $self->display_node;
  my $sub = $self->Sub_Section;
  while( $sub->File !~ /NULL/ ) {
    $sub->display_node;
    $sub = $sub->Sub_Section;
  }
  print "\n\nEND DISPLAY OF NODES:\n";
}

sub display_titles {
#
# Display the titles of the tree
#
  my $self = shift;

  my $title = $self->Title;
  print "$title\n";
  my $sub = $self->Sub_Section;
  while( $sub->File !~ /NULL/ ) {
    my $title = $sub->Title;
    print "$title\n";
    $sub = $sub->Sub_Section;
  }
}
sub Get_Sub_Section_Links( )
#
# Get the links to the main sub-sections from this section
#
{
  my $self = shift;
  my $section = shift;

  my $sect_lev = $section->Level;
  my $sub = undef;
  $sub = $section->Sub_Section;
  #
  # If sub-section defined
  #
  if ( defined($sub) && ($sub->File !~ /NULL/) ) {
    my $level = $sub->Level;
    #
    # If next section is a sub-section of this section
    #
    if ( $level > $sect_lev ) {
      my $sub_level = $level;
      #
      # Loop over all sections of this level until we reach a section that
      # is of the same level as the main section.
      #
      my $sect = $self;
      my $sub_level_one = $sub_level;
      while( (defined($sub)) && ($sub_level > $sect_lev) ) {
        #
        # Save each section that is of this level
        #
        $level = $sub->Level;
        if ( $level == $sub_level_one ) {
          my $title = $sub->Title;
          my $file = $sub->File;
          $sect->{TITLE} = $title;
          $sect->{FILE} = $file;
          $sect->{LEVEL} = 2;
          $sect->{SUB_SECTION} = Section->new( );
          $sect = $sect->{SUB_SECTION};
        }
        $sub = $sub->Sub_Section;
        $sub_level = $sub->Level;
      }
    }
  }
}

sub Make_Section_Link
#
# Figure out the section links based on the history of the up refs
# And write this section history to the given file handle.
#
{
  my $self = shift;
  my $fh   = shift;
  my $title= shift;
  my $file = shift;
  my $up   = shift;
  my $sub = undef;
  #
  # If first file get the original title as the first link
  #
  if ( $self->File =~ /NULL/ ) {
    $self->Fill_Section( $title, $file, $up );
  }
  else {
    #
    # If a new section redo the section links
    #
    if ( ! defined( $up ) ) {
      $self->Fill_Section( $title, $file, $up );
    }
    #
    # If next file
    #
    elsif ( $self->Sub_Section->File =~ /NULL/ ) {
      $self->Fill_Section( $title, $file, $up );
    }
    #
    # If a sub-section of a previous section
    #
    else {
      #
      # First find what level the new section is so we can find it's place
      #
      $sub = $self->Sub_Section;
      my $level = Title_Level( $title );
      while( defined($sub->Sub_Section) && 
      ($sub->Sub_Section->File !~ /NULL/) && 
      ($sub->Sub_Section->{LEVEL} < $level) ) {
        $sub = $sub->Sub_Section;
      }
      $level = $sub->Level;
      $sub->Fill_Section( $title, $file, $up );
    }
  }
  my $toc = "False";
  $self->write_HTML( $fh, $toc  );
}

sub Make_TOC_Link
#
# Figure out the table-of-contents links based on the history of the up refs
#
{
  my $self = shift;
  my $title= shift;
  my $file = shift;
  my $up   = shift;

  my $sub = undef;
  #
  # If first file get the original title as the first link
  #
  if ( $self->File =~ /NULL/ ) {
    $self->Fill_Section( $title, $file, $up );
  }
  #
  # If next file
  #
  elsif ( $self->Sub_Section->File =~ /NULL/ ) {
    $self->Fill_Section( $title, $file, $up );
  }
  #
  # If below first or second link
  #
  else {
    $sub = $self->Sub_Section;
    while( defined($sub->Sub_Section) && ($sub->Sub_Section->File !~ /NULL/) ) {
      $sub = $sub->Sub_Section;
    }
    $sub->Fill_Section( $title, $file, $up );
  }
}


1  # To make use or require happy
