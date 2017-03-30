#!/usr/bin/env perl
#
#       html2word.pl            Erik Kluzek             Oct/7/97
#                               NCAR/CGD/CSM
#                               erik@ncar.ucar.edu
#                               (303)497-1326
#
#	Convert the CAM users guide HTML source files to a file that can be
#	read in Word for a paper-based copy.
#
#
use strict;
use diagnostics;
use lib ( ".", "perl_mod" );
use CAM_Web_File;
use html2w;

my $outname = "cam_word";
my $outnum = 0;
my $word = HTML2Word->new( "index.shtml" );
$word->open_out( $outname, $outnum );

#
# Loop over the files in the correct order...
#
my $next = $word->File;
my $nfile = 0;
while ( defined( $next ) ) {
  my $file = $next;
  $nfile++;
  print "File: ($nfile) $file\n";
  $word->strip_file_headers( $file );
  $next = $word->Next;
}
#
# Realign (shorten lines, get the comments lined up etc.) the output file
#
$word->close_out;
