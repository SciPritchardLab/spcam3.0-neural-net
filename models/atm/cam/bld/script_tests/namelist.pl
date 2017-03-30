#!/usr/bin/env perl
#
#	namelist.pl			Erik Kluzek
#
#	Unit test of the namelist base-class.
#
#	$Id: namelist.pl,v 1.1.6.2 2002/07/31 20:00:16 erik Exp $
#

use strict;
#use diagnostics;
use lib "..", "../../../bld";

use namelist;

%main::CCMEXP = {};
%main::LSMEXP = {};
$main::CCMEXP{'NCDATA'}  = "\'/CCM/csm/input/atm/ccm3/ccm3/SEP1.T42L26.112000.nc\'";
$main::CCMEXP{'nElapse'} = -1;
my $file = "nl.initial";
my $nl = namelist->new( "CCMEXP", $file, \%main::CCMEXP );
print "Here is the namelist\n";
$nl->print;
print "convert_case\n";
$nl->convert_case;
print "print again\n";
$nl->print;
print "Here is what the keys to %main::CCMEXP looks like now\n";
my @list = keys(%main::CCMEXP);
print "@list\n";
print "change\n";
$nl->change;
print "check\n";
$nl->checkstring( "ncdata" );
print "Write\n";
$nl->Write( "." );
system( "cat $file" );
print "Append a new namelist on the end\n";
my $lmnl = namelist->new( "LSMEXP", $file, \%main::LSMEXP );
$main::LSMEXP{'finiDAT'}  = "\'/CCM/csm/input/lnd/lsm1/CCM3.10LSMICSEP1.072000.nc\'";
$lmnl->convert_case;
$lmnl->Write( ".", "Append" );
system( "cat $file" );
print "Now parse the namelist\n";
%main::CCMEXP = {};
%main::LSMEXP = {};
$nl->parse( $file );
$nl->print;
$lmnl->parse( $file );
$lmnl->print;
print "Now parse a more complex namelist:\n";
%main::CCMEXP = {};
%main::LSMEXP = {};
open( NAMELIST, ">$file" ) || die "ERROR: can not open $file\n";
print NAMELIST <<EOF;
&ccmexp
 string = 4*"4times"
 string2 = "1times", 2*"2times","stuff", 		
	3*"3times"
 fourones = 4*1
 two5sa4a3andtwo8s = 2*5, 4, 3, 2*8
 one1two5s = 1,2*5
 a = 'b'
 absems_data    = '/fs/cgd/csm/inputdata/atm/ccm3/abs_ems_factors_fastvx.052001.nc'
 ncdAta = "as'dfasdfasdf'asdf", STUFF = -1, !comment
 adiabatic      = .false.,
 bndtvg         = '/fs/cgd/csm/inputdata/atm/ccm3/noaamisc.r8.nc'
 BNDTVo         = '/fs/cgd/csm/inputdata/atm/ccm3/noaao3.1990.21999.nc'
 bNDTvs         = '/fs/cgd/csm/inputdata/atm/ccm3/T21M5079.nc'
 caseid         = 'testaixeulT21L18'
 cTITle         = 'testaixeulT21L18 &;:[]{}+!=\@:# $|\^(*)~dyn: eul res: T21' !comment
 dif2           = 250000, dIF4           = 2e+16, fincl1(1)      = 'T:I', fincl1(2)      = 'PS:I'
 phys_grid%n    = 1, phys_grid%thing = 1e+16_8, pp_grid%value = -1.123456E-349_dbl_kind,
 fincl2(1)      = 'T:A'	!comment
 fincl2(10)     = 'Z300:A'
 fincl2(11)     = 'Z500:A'
 fincl2(12)     = 'Z700:A'
 fincl2(2)      = 'PS:A'
 fincl2(3)      = 'PSL:A'
 fincl2(4)      = 'U200:A'
 fincl2(5)      = 'U850:A'
 fincl2(6)      = 'V200:A'
 fincl2(7)      = 'V850:A'
 fincl2(8)      = 'T300:A'
 fincl2(9)      = 'T850:A'
 fexcl1         = 'T', 'Q', 'U',
                  'Z', 'W', 'L',
                  'N', 'Z', 'W', 'Y', 'L', 'M', 'N'
                  'O', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', 'ZZ'
 fexcl2         = 'T' 'Q' 'U'


     			
 fexcl6         = 'TREFAT:I', 'TREFHA:I', 'aREFHT:I' 'TaEFHT:I'
                  'TREFBT:I', 'TREFHB:I', 'bREFHT:I' 'TbEFHT:I', 'PSL:X'
                  'TREFCT:I', 'TREFHc:I','cREFHT:I' 'TcEFHT:I',
                  'TREFDT:I', 'TREFHd:I', 'dREFHT:I' 'TdEFHT:I'
                  'TREFET:I', 'TREFHe:I', 'eREFHT:I' 'TeEFHT:I'
                  'TREFFT:I', 'TREFHf:I', 'fREFHT:I' 'TfEFHT:I'
                  'TREFGT:I', 'TREFHg:I','gREFHT:I' 'TgEFHT:I', 
                  ,  'TREFHT:I', 'TREFHh:I', 'hREFHT:I' 'ThEFHT:I',
                  ,'TREFIT:I', 'TREFHi:I', 'iREFHT:I' 'TiEFHT:I',
                  ,"TREFJT:I", 'TREFHj:I', "jREFHT:I" 'TjEFHT:I'
                  'TREFKT:I', 'TREFHk:I','kREFHT:I' 'TkEFHT:I'
                  'TREFLT:I', 'TREFHl:I', 'lREFHT:I' 'TlEFHT:A', 'OBLIQ_PSL:I',
                  'SPECHUM_1000:I', 'RELHUM:A', 'PSL_1000:A', 'U_850:I'
                  'OMEGA_850MB:X,q', '5MEGA_850MB:X' 'SPECHUM_1001:I'
                  '1MEGA_850MB:X,q', '4MEGA_850MB:X' 'SPECHUM_1002:I'
                  '2MEGA_850MB:X,q','3MEGA_850MB:X' 'SPECHUM_1003:I'
                  "3MEGA_850MB:X", '2MEGA_850MB:X' 'SPECHUM_1004:I'
                  "4MEGA_850MB:X", '1MEGA_850MB:X' 'SPECHUM_1005:I'
                  '5MEGA_850MB:X', '0MEGA_850MB:X' 'SPECHUM_1006:I'
                  '6MEGA_850MB:X','9MEGA_850MB:X' 'SPECHUM_1007:I'
                  '7MEGA_850MB:X', '8MEGA_850MB:X' 'SPECHUM_1008:I'
 thingy1        = 1900,1901, 1902, 1903 1904, 1905, 1906, 1907, 1908, 1909
                  1911,1911, 1912, 1913 1914, 1915, 1916, 1917, 1918, 1919
                  1921,1921, 1922, 1923 1924, 1925, 1926, 1927, 1928, 1929
                  1931,1931, 1932, 1933 1934, 1935, 1936, 1937, 1938, 1939
                  1941,1941, 1942, 1943 1944, 1945, 1946, 1947, 1948, 1949
                  1951,1951, 1952, 1953 1954, 1955, 1956, 1957, 1958, 1959
                  1961,1961, 1962, 1963 1964, 1965, 1966, 1967, 1968, 1969
                  1971,1971, 1972, 1973 1974, 1975, 1976, 1977, 1978, 1979
                  1981,1981, 1982, 1983 1984, 1985, 1986, 1987, 1988, 1989
                  1991,1991, 1992, 1993 1994, 1995, 1996, 1997, 1998, 1999
 logical = 4*.false.,2*.true.
 logical2 = .true.,2*.false.
 real =4*5.6_8,2*0.0
 real2 = +9.999D00,4*5.6,2*0.0
 complex = 100*(4,5.0e-10)
 complex2 = (-5.2345d+23,+5.98765E-100),2*(4,5.0e-10)
 irt            = 0
 iyear_ad       = 1950
 ndens(1)       = 1
 ndens(2)       = 1
 nestep         = 30
 nhtfrq(1)      = 1
 nhtfrq(2)      = 10
 nnbdat         = 1231
 nsrest         = 0
 readtrace      = .false.
 trace_gas      = .true.
 hist_mfilt     =  1,1
/
EOF
close( NAMELIST );
$nl->parse( $file );
print "Print the complex namelist out\n";
$nl->print;
