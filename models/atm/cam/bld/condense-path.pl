#!/usr/bin/env perl                                                   
#
#	condense_path.pl			Brian Eaton
#
# 	Remove "dir/.." from a directory path.                              
#
#	$Id: condense-path.pl,v 1.1.2.2 2002/02/09 01:19:01 eaton Exp $
#
use Env qw(MODEL_CFGDIR);
use strict;
use lib (".","$MODEL_CFGDIR");   # List of where to look for the Perl modules
use CAM;

my $path = shift @ARGV;                                                  
my $newpath = &CAM::condense_path( $path );
print $newpath;
