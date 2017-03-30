#!/usr/bin/env perl
#
#       word2html.pl           Erik Kluzek             
#                              NCAR/CGD/CSM
#                              erik@ncar.ucar.edu
#                              (303)497-1326
#
# Convert Word-html paper-copy version back to web-based version.
#
# $Id: word2html.pl,v 1.18.2.3 2002/12/04 16:18:32 erik Exp $
#
use strict;
use diagnostics;
use lib ( ".", "perl_mod" );
use html2w;
use CAM_w2html;

my $html = CAM_Word2HTML->new( );
my $file;
#
# Convert the table of contents
#
#my @pre   = ( "toc" );
#my @files = ( "table_of_contents.html" );
#my @outfile = ( undef );
#my $file = undef;
#my $i = undef;
#for( $i = 0;  $i <= $#pre; $i++ ) {
#  my $word = HTML2Word->new( $files[$i] );
#  $word->open_out( $pre[$i], 0 );
#  $outfile[$i] = $word->Outfile;
#  $file = $word->File;
#  print "Convert $file\n";
#  $word->strip_file_headers( $file );
#  $word->close_out;
#}
#
# Loop over all cam_word.html files
#
my @first_file = ( "index.shtml" );
while( glob("cam_word0.html") ) {
  my $infile = $_;
  print "infile: $infile\n";
  my $file = shift @first_file;
  $html->open_in( $infile, $file );
  $html->put_back_file_headers;
}
print "Now go through and put the sub-section headers in the files\n";
my $web  = $html->Web;
my $sect = $web->Table_Of_Cont;
my $prev = $sect->File;
$sect = $sect->Sub_Section;
while( defined($sect) && ($sect->File !~ /NULL/) ) {
  $file = $sect->File_Name;
  my $sublink = Section->new;
  $sublink->Get_Sub_Section_Links( $sect );
  my $level = $sublink->Level;
  if ( ($file =~ /table_of_contents/) || ((defined($level)) && ($file !~ /$prev/)) ) {
    print "Add sub-links to file: $file\n";
    $web->link( $sublink );
    my $word = HTML2Word->new( $file );
    $word->open_out( "cam_tmp", 0 );
    my $out_file = $word->File;
    $word->strip_file_headers( $out_file );
    $word->close_out;
    my $infile = $word->Outfile;
    $html->open_in( $infile, $prev );
    $html->put_back_file_headers;
  }
  $sublink->DESTROY;
  $prev = $file;
  $sect = $sect->Sub_Section;
}
