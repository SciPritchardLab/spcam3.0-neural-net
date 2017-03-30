#!/usr/bin/env perl
#
#	realign.pl				
#      Erik Kluzek
#	Feb/28/98
#
#	Realign the HTML comments from the given input HTML file.
#	This is to ensure that any other operations that look for HTML
#	comments will work properly.  In some cases those operations may
#	not function correctly if the comments are not aligned correctly.
#
#	Also shortens long lines so that they fit within a reasonable line length
#	near 120 characters.
#
#	$Id: realign.pl,v 1.4.2.2 2002/05/30 17:40:36 erik Exp $
#
use strict;

sub truncate
{
#
# truncate( string )
#
  my $line = shift;
  my $new_line = "";
  my $long_line_len = 120;

  if ( (length( $line ) > $long_line_len) && ($line !~ /<TITLE>/) && ($line !~ /<!/) && ($line !~ /<H[123456789]/)) {
    my $begin = $long_line_len - 20;
    my $end = $long_line_len + 20;
    my $string = $line;
    if ( $string =~ /^(.{$begin,$end})([ \t]+.*$)/ ) {
      $new_line = $1;
      $string = $2;
      while( $string =~ /(^.{$begin,$end})([ \t]+.*$)/ ) {
        $new_line = "$new_line\n$1";
        $string = $2;
      }
      $new_line = "$new_line\n$string\n";
    }
    else {
      $new_line = $line;
    }
  }
  else {
    my $string = $line;
    if ( $string =~ /(<!.*)/ ) {
      chop( $string ); # remove new-line
      my @html =  split /</, $line;
      my $i = undef;
      $new_line = "";
      foreach $i (@html) {
        if ( $i =~ /.+/ ) {
          $new_line = "$new_line<$i";
          if ( $i != $html[$#html] ) {
            $new_line = "$new_line\n";
          }
        }
      }
    }
    elsif ( $string =~ /(.+)(<TITLE>*)/ ) {
      $new_line = "$1\n$2\n";
    }
    elsif ( $string =~ /(.+)(<H[123456789]>*)/ ) {
      $new_line = "$1\n$2\n";
    }
    else {
      $new_line = $string;
    }
  }
  return( $new_line );
}

  my $file = shift;
  if ( ! defined($file) ) {
    die "realign file not given\n";
  }
  open( FILE, "<$file" ) || die "Could not open file $file";
  open( OUT, ">$file.tmp" ) || die "Could not open file $file";
  while( <FILE> ) {
    $_ = &truncate( $_ );
    print OUT $_;
  }
  close( FILE );
  close( OUT );
  system( "/bin/mv $file.tmp $file" );
