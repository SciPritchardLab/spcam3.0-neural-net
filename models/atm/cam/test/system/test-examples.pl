#!/usr/bin/env perl
#
# Test script to run the Users-Guide examples on the NCAR IBM (blackforest)
#
#-----------------------------------------------------------------------
# Batch options for machine with loadleveler (IBM SP) (blackforest)
# submit with
#        llsubmit test-examples.pl
# On IBM with batch submission the number of nodes won't change
# even when the PE configuration is changed.
#-----------------------------------------------------------------------
# Name of the que (CHANGE THIS if needed)
# @ class       = csl_pr
# Number of nodes (CHANGE THIS if needed)
# @ node        = 8
# Switch to use (CHANGE THIS if needed)
# @ network.MPI = csss,not_shared,us
# @ output      = testex.aix.log
# @ error       = testex.aix.err
# @ node_usage  = not_shared
# @ job_type    = parallel
# @ tasks_per_node = 1
# Export all Environment variables
# @ environment = COPY_ALL
# @ queue
#-----------------------------------------------------------------------
use Env qw(CSMDATA LOGNAME);
use Cwd;

#
# Specifics for running at a given site
#
my $resume = shift(@ARGV);
my $SPMDRUN_CMND = "poe";
my $CASE_DIR = "/ptmp/$LOGNAME";
my $GNUMAKE = "gmake";
my $first_test = 1;
my $ncdump = "ncdump -v date,datesec";

if ( $CSMDATA !~ /.+/ ) { $CSMDATA = "/fs/cgd/csm/inputdata"; }

#
# Setup the namelists options that will be used for each example
#
my $pwd = cwd();
my $configure = "$pwd/../../bld/configure";
my $build_namelist = "$pwd/../../bld/build-namelist";
my $file4 = "namelist_input_run4";
my $file5 = "namelist_input_run5";
my $file6 = "namelist_input_run6";
my $file7 = "namelist_input_run7";
my $file8 = "namelist_input_run8";
my $config = "config_cache.xml";
my %options = ( test1=>"-config obj/$config",
                test2=>"-config obj/$config -case run02  -namelist  \"&camexp start_ymd=010901 start_tod=0 stop_ymd=010902 stop_tod=0 mss_irt=0 /\"",
                test3=>"-config obj/$config -case run02 -runtype restart  -namelist '&camexp nelapse=-30 /'", 
                test4=>"-config obj/$config -case run04 -infile $file4",
                test5=>"-config obj/$config -case run05 -runtype branch -infile $file5",
                test6=>"-config obj/$config -case run06 -infile $file6",
                test7=>"-config obj/$config -case run07 -infile $file7",
                test8=>"-config obj/$config -case run08 -infile $file8",
                test9=>"-config obj/$config -case run09 -namelist \"&camexp ncdata='$CSMDATA/atm/cam2/inic/gaus/JAN1.T42L26.navytopo.FRAC.ross.c020305.nc', start_ymd=101 / &clmexp finidat=' ' /\""
              );
#
# Set up the directories to use
#
my $run = "$CASE_DIR/test_examples";
if ( ! -d $run ) { mkdir( $run, 0755 ) || die "ERROR:: Can not create $run directory\n"; }
my $build = "$run/obj";
if ( ! -d $build ) { mkdir( $build, 0755 ) || die "ERROR:: Can not create $run directory\n"; }
print "Change to $build directory\n";
chdir( $build ) || die "ERROR: Can not change to $build directory\n";
#
# Figure out up-case name for MSS
#
my $MSS = $LOGNAME;
$MSS =~ tr/a-z/A-Z/;
#
# Configure and build
#
my $conf = "-cam_bld $build -cam_exedir $run -cache $config";
if ( ! -f "misc.h" ) {
  print "$configure $conf\n";
  system( "$configure $conf" );
  if ( $? != 0 ) {
    die "ERROR in configure\n\n";
  }
}
print "$GNUMAKE > compile.log  2>&1\n";
system( "$GNUMAKE > compile.log  2>&1" );
if ( $? != 0 ) {
  die "ERROR in gmake\n\n";
}
print "Change to $run directory\n";
chdir( $run ) || die "ERROR: Can not change to $run directory\n";
#
# Run through all tests at default resolution
#
foreach my $test ( sort(keys(%options)) ) {
  if ( defined($resume) ) {
    if ( $test lt $resume ) { print "Skip test: $test\n"; next; }
  }
  print "Test: $test\n";
  if (      $test eq "test4" ) {
    open( FILE, ">$file4" ) || die "ERROR:: Could not open $file4\n";
    print FILE << "EOF";
&camexp 
 ctitle = 'Paleo-climate run: 4050 BC' 
 eccen = 0.01868126 
 obliq = 24.10538 
 mvelp = 0.87 
 co2vmr = 2.80e-4 
 ch4vmr = 0.700e-6 
 n2ovmr = 0.275e-6 
 f11vmr = 0. 
 f12vmr = 0. 
 scon = 1.363e6 
 / 
EOF
    close(FILE);
  } elsif (      $test eq "test5" ) {
    open( FILE, ">$file5" ) || die "ERROR:: Could not open $file5\n";
    print FILE << "EOF";
&camexp 
 caseid = 'run05' 
 nrevsn = '/$MSS/csm/camrun/atm/rest/camrun.cam2.r.0000-09-02-00000' 
 nelapse = -30 
 nhtfrq = 72,18,18 
 ndens = 2,2,2 
 mfilt = 10,30,30 
 fincl2 = 'PRECL:I','PRECC:I','CLOUD:I' 
 fincl3 = 'PRECL:M','PRECC:M','CLOUD:M' 
 fexcl1 = 'PRECL','PRECC','CLOUD' 
 mss_wpass = 'mypass' 
 / 
 &clmexp 
 nrevsn = '/$MSS/csm/camrun/lnd/rest/camrun.clm2.r.0000-09-02-00000' 
 / 
EOF
    close(FILE);
  } elsif ( $test eq "test6" ) {
    open( FILE, ">$file6" ) || die "ERROR:: Could not open $file6\n";
    print FILE << "EOF";
 &camexp 
 mss_irt = 1825 
 mfilt = 10 
 nelapse = -20 
 nhtfrq = 0,36 
 mfilt = 1,2 
 fincl1 = 'PRECSL','PRECSC' 
 fincl2 = 'T','U' 
 fexcl2 = 'T850' 
 hfilename_spec(2) = 'h%t.%y-%m-%d.nc' 
 avgflag_pertape(2) = 'X' 
 / 
EOF
    close(FILE);
  } elsif ( $test eq "test7" ) {
    open( FILE, ">$file7" ) || die "ERROR:: Could not open $file7\n";
    print FILE << "EOF";
&camexp 
 ctitle = 'multiyear SST (AMIP)dataset' 
 bndtvs = '$CSMDATA/atm/cam2/sst/sst_HadOIBl_bc_64x128_1949_2001_c020411.nc'
 start_ymd = 19790901 
 start_tod = 0 
 nelapse = -180
 sstcyc = .false. 
 / 
EOF
    close(FILE);
  } elsif ( $test eq "test8" ) {
    open( FILE, ">$file8" ) || die "ERROR:: Could not open $file8\n";
    print FILE << "EOF";
&camexp 
 nrefrq = 0 
 mss_irt = 0 
 linebuf = .TRUE. 
 ndens = 1 
 mfilt = 10 
 nhtfrq = 1 
 nelapse= 10 
 / 
EOF
    close(FILE);
  }
  if ( ($test cmp "test$first_test") >= 0 ) {
    print "Namelist option: $options{$test}\n";
    print "$build_namelist -o $test.nl $options{$test}\n";
    system( "$build_namelist -o $test.nl $options{$test}" );
    if ( $? != 0 ) {
      die "ERROR in build-namelist $test\n\n";
    }
    print "$SPMDRUN_CMND cam < $test.nl > $test.log\n";
    system( "$SPMDRUN_CMND cam < $test.nl > $test.log" );
    if ( $? != 0 ) {
      die "ERROR in cam $test\n\n";
    }
    #
    # Validate and verify output
    #
    if ( $test eq "test5" ) {
      my ($file0) = glob("run05.cam2.h0.*.nc");
      if ( ! -r $file0 ) { die "ERROR:: Can not read in $file0\n"; }
      my ($file1) = glob("run05.cam2.h1.*.nc");
      if ( ! -r $file1 ) { die "ERROR:: Can not read in $file1\n"; }
      my ($file2) = glob("run05.cam2.h2.*.nc");
      if ( ! -r $file2 ) { die "ERROR:: Can not read in $file2\n"; }
      `$ncdump $file0 > $file0.dump`;
      if ( $? != 0 ) {
        die "ERROR:: Problem doing ncdump on $file0\n";
      }
      `$ncdump $file1 > $file1.dump`;
      if ( $? != 0 ) {
        die "ERROR:: Problem doing ncdump on $file1\n";
      }
      `$ncdump $file2 > $file2.dump`;
      if ( $? != 0 ) {
        die "ERROR:: Problem doing ncdump on $file2\n";
      }
      &validate( "time = UNLIMITED ; .. .(\d+) currently", "10", "$file0.dump" ) ;
      &validate( "time = UNLIMITED ; .. .(\d+) currently", "30", "$file1.dump" ) ;
      &validate( "time = UNLIMITED ; .. .(\d+) currently", "30", "$file2.dump" ) ;
      &verify( "mswrite  -t  365 -w mypass", "$test.log" );
      &verify( " date = 923, 924", "$file0.dump" );
      &verify( " datesec = 0, 0",  "$file0.dump" );
      &verify( " date = 924, 925", "$file1.dump" );
      &verify( " datesec = 64800, 0",  "$file1.dump" );
      &verify( " date = 924, 925", "$file2.dump" );
      &verify( " datesec = 64800, 0",  "$file2.dump" );
      foreach my $var ( ("PRECL", "PRECC", "CLOUD") ) {
        &verify_none( "float $var", "$file0.dump" );
        &verify( "$var:cell_method = .time: minimum.", "$file2.dump" );
        &verify_none( "$var:cell_method", "$file1.dump" );
        &verify( "float $var", "$file1.dump" );
        &verify( "float $var", "$file2.dump" );
      }
    }
    if ( $test eq "test7" ) {
      # Validate that start date is expected value
      &validate( "Start date\s+([^ ]+)", "19790901", "$test.log" );
      # Verify that January of 1980 was read in for SST data
      &verify( "Read sst for date (yyyymmdd)  19800116", "$test.log" );
    }
    if ( $test eq "test6" ) {
      my ($file1) = glob("h?.????-??-??.nc");
      if ( ! -r $file1 ) { die "ERROR:: Can not read in $file1\n"; }
      `$ncdump $file1 > $file1.dump`;
      if ( $? != 0 ) {
        die "ERROR:: Problem doing ncdump on $file1\n";
      }
      &verify( "[A-Z0-9_-]+:cell_method = .time: maximum.", "$file1.dump" );
      &verify( "mswrite  -t 1825", "$test.log" );
    }
    if ( $test eq "test8" ) {
      my ($file0) = glob("run08.cam2.h0.*.nc");
      if ( ! -r $file0 ) { die "ERROR:: Can not read in $file0\n"; }
      `$ncdump $file0 > $file0.dump`;
      if ( $? != 0 ) {
        die "ERROR:: Problem doing ncdump on $file0\n";
      }
      &verify_none( "float ", "$file0.dump" );
      &verify( "double [A-Z0-9_-]+", "$file0.dump" );
      &verify_none( "mswrite", "$test.log" );
      &verify_none( "WRITE_REST_PFILE", "$test.log" );
    }
  }
}
#-----------------------------------------------------------------------
# Now run test2 at 32x64_T21 resolution
#-----------------------------------------------------------------------
#
# Set up the directories to use
#
$run = "$run.32x64";
if ( ! -d $run ) { mkdir( $run, 0755 ) || die "ERROR:: Can not create $run directory\n"; }
my $build = "$run/obj";
if ( ! -d $build ) { mkdir( $build, 0755 ) || die "ERROR:: Can not create $run directory\n"; }
print "Change to $build directory\n";
chdir( $build ) || die "ERROR: Can not change to $build directory\n";
#
# Configure and build
#
print "Change to $build directory\n";
chdir( $build ) || die "ERROR: Can not change to $build directory\n";
my $conf = "-cam_bld $build -cam_exedir $run -cache $config -res 32x64";
print "$configure $conf\n";
system( "$configure $conf" );
if ( $? != 0 ) {
  die "ERROR in configure\n\n";
}
print "$GNUMAKE > compile.log  2>&1\n";
system( "$GNUMAKE > compile.log  2>&1" );
if ( $? != 0 ) {
  die "ERROR in gmake\n\n";
}
print "Change to $run directory\n";
chdir( $run ) || die "ERROR: Can not change to $run directory\n";
#
# Build namelist and run
#
my $test = "test2";
print "$build_namelist $options{$test}\n";
system( "$build_namelist $options{$test}" );
if ( $? != 0 ) {
  die "ERROR in build-namelist $test\n\n";
}
print "$SPMDRUN_CMND cam < namelist > $test.log\n";
system( "$SPMDRUN_CMND cam < namelist > $test.log" );
if ( $? != 0 ) {
  die "ERROR in cam $test\n\n";
}

sub verify {
#
# Verify that a given string exists in the input file
#
  my $pattern = shift;
  my $filename = shift;

  open( FILE, "<$filename" ) || die "ERROR:: Can not open file: $filename\n";
  my $match = 0;
  while( $_ = <FILE> ) {
    if ( /$pattern/ ) { $match = 1; }
  }
  close( FILE );
  if ( ! $match ) {
    die "ERROR:: Could not find $pattern in $filename\n";
  }
}

sub verify_none {
#
# Verify that a given string does NOT exist in the input file
#
  my $pattern = shift;
  my $filename = shift;

  open( FILE, "<$filename" ) || die "ERROR:: Can not open file: $filename\n";
  my $match = 0;
  while( $_ = <FILE> ) {
    if ( /$pattern/ ) { $match = 1; }
  }
  close( FILE );
  if ( $match ) {
    die "ERROR:: Found $pattern in $filename when it should not be there\n";
  }
}

sub validate {
#
# Validate that a given pattern exists in the input file and the value
# of a variable matches the expected value.
#
  my $pattern = shift;
  my $expect = shift;
  my $filename = shift;

  open( FILE, "<$filename" ) || die "ERROR:: Can not open file: $filename\n";
  my $match = 0;
  while( $_ = <FILE> ) {
    if ( /$pattern/ ) { 
      if ( $1 ne "$expect" ) {
        die "ERROR:: Found $pattern in $filename but does not match expected value of $match (does equal $1)\n";
      }
    }
  }
  close( FILE );
  if ( $match ) {
    die "ERROR:: Found $pattern in $filename when it should not be there\n";
  }
}
