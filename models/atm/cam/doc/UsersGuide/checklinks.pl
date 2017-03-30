#!/usr/bin/env perl
#
#       checklinks.pl           Erik Kluzek             Apr/14/98
#                               NCAR/CGD/CSM
#                               erik@ncar.ucar.edu
#                               (303)497-1326
#
#
# Go through the cam_word0.html file and check all the HTML link
# references.  Ensure that the files and targets referenced actually
# exist.  If not print an error, and also output the error and the list
# of links to a set of output files (cam_doc*.out).
#
#	Version: $Id: checklinks.pl,v 1.4.2.1 2002/03/11 21:41:40 erik Exp $
#
use strict;
use lib ( ".", "perl_mod");
use Check_HTML_Link;

my $file = "cam_word0.html";
open( FILE, "<$file" ) || die "could not open file: $file";

my $clink = Check_HTML_Link->new( $file );
print "Loop over $file\n";
while( <FILE> ) {
  $clink->check( $_ );
} 
print "Finished file: \n"
