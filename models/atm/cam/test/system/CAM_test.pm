#
#	CAM_test.pm			Erik Kluzek
#
#	Perl5 Object module extending the CAM_run.pm module
#	to handle testing of the Community Atmospheric Model (CAM).
#
#	Methods:
#
#	new ---------------- Constructor
#	process_args ------- Process main script command line arguments.
#	setup_tests -------- Setup the tests that will be run.
#	list_tests --------- List the tests performed by the test suite.
#	usage -------------- Terminate and return proper syntax for arguments.
#	setup_directories -- Setup directories.
#	remote_baseline_errg The filename of the remote error-growth file.
#	dynamics ----------- List of dynamics to run tests on.
#	tests -------------- List of tests to run.
#	config ------------- Return configuration setting for given test.
#	exec --------------- Extend CAM_run exec method so that commands aren't
#				executed if you are continuing from a previous run and
#				haven't caught up to the point where you started.
#	build_namelist ----- Extend CAM_run build_namelist method so that it's only
#				executed if you aren't continuing from a previous run.
#	do_compare --------- Return true if you want to compare to another program library.
#	do_continue -------- Return true if the continue option is set.
#	do_nofail ---------- Return true is the nofail option is set.
#	do_remote_errg ----- Return true if you want error-growth tests to compare to 
#				remote baseline.
#	numerically -------- Use with "sort" to sort values numerically.
#	completed ---------- Return true if this part of the script was already completed.
#	compare_files ------ Compare history files of two simulations (uses cprnc)
#	report_compare ----- Summary report on success of file comparisions.
#	mark --------------- Mark to describe where at in process.
#	start_log ---------- Initialize the test-log file.
#	read_mark ---------- Read the continuation mark.
#	write_mark --------- Write the continuation mark.
#	continuation_mark -- Mark how far the script has gone.
#	lower_PEs ---------- Lower the number of CPU's being used (as a test).
#	reset_PEs ---------- Reset the number of CPU's being used.
#	die ---------------- Terminate script or continue if no-fail option set.
#	version ------------ Get the version from the ChangeLog.
#	save_timing -------- Save timing information for a test.
#
#	$Id: CAM_test.pm,v 1.6.2.18 2003/06/16 14:55:52 erik Exp $
#
use 5.004;   # Use at least this version of perl
use strict;
#use diagnostics;
use Cwd;

package CAM_test;
@CAM_test::ISA = "CAM_run";
use CAM_run;
use cam_timing;
#
# Use all of the CAM_run methods and data
# Override the constructor, process_args, and usage methods.
#

sub new {
#
# Constructor
#
  my $class = shift;

  my $self = $class->SUPER::new;
  my $OS = $self->{'os'};
  $self->{'clean'} = "yes";           # By default clean the old directories out
  $self->{'compare'} = "no";          # Compare to control source
  $self->{'rootdir'}  = undef;        # Original CAMROOT if it changes
  $self->{'cfgdir'}  = undef;         # Original MODEL_CFGDIR if it changes
  $self->{'cam_cont_dir'} = undef;    # Directory with CAM model source to compare
  $self->{'continue'} = "no";         # If you want to continue a previous script run
  $self->{'completed'} = "yes";       # If this part has already been completed or not
  $self->{'nofail'} = "no";           # Don't die on errors, write log and continue
  $self->{'errors'} = 0;              # Number of errors when running in nofail mode
  $self->{'status'} = "ran";          # Status of test (ran or FAIL)
  $self->{'continue_file'} = undef;   # filename that keeps track of how far script has gone
  $self->{'testlog'} = undef;         # Log of results
  $self->{'firsttest'} = undef;       # First test to run
  $self->{'lasttest'} = undef;        # Last test to run
  $self->{'begin'} = "yes";           # Mark that script has just begun
  $self->{'remote_machine'}  = undef; # Remote machine to use for error-growth
  $self->{'remote_platform'} = undef; # Remote platform to use for error-growth
  $self->{'remote_lab'} = undef;      # Remote lab to use for error-growth
  $self->{'dynamics'} = "all";        # Dynamics to run tests on
  $self->{'run_res'} = {};            # Hash of run resolutions
  $self->{'run_lev'} = {};            # Hash of run vertical levels
  $self->{'debug_config'} = {};       # Hash of configuration settings for debug tests
  $self->{'runtypes'} = {};           # Hash of run-types for restart tests
  $self->{'nestep'} = {};             # Hash of nestep for restart tests
  $self->{'roundoff'} = undef;        # Value to use for roundoff tests
  $self->{'errgrotests'} = undef;     # Error growth tests
  $self->{'physics_tests'} = undef;   # physics test descriptions
  $self->{'control_tests'} = undef;   # control test descriptions
  $self->{'timing'} = {};             # Hash of timing information
  $self->{'twod_yz'} = "FALSE";       # 2D decomposition for FV dynamics
  $self->{'use_unique_id'} = "no";    # unique job ID for concurrent runs of 
  $self->{'unique_id'} = undef;       # test-model on the same machine
  #
  # Env variables
  #
  my @more_env_list = (
      "EXEDIR",          # Executable directory name
      "OCEANMODEL",      # Ocean model to use (dom or som)
      "PHYSICS",         # Physics to use (cam1 or ccm366)
      "MODEL_CPRDIR",    # cprnc directory to compare history files
      "CONT_ROOTDIR"     # Control root directory
  );
  my $env_ref  = $self->{'env_list'};
  my @env_list = @$env_ref;
  foreach my $env ( @more_env_list ) {
    $self->{$env} = undef;
  }
  push( @env_list, @more_env_list );
  $self->{'env_list'} = \@env_list;
  bless( $self, $class );
  $self->import_env;
  # default file paths
  my $SCRIPT_DIR = $self->env( "SCRIPT_DIR" ); 
  $self->{'testlog'} = "$SCRIPT_DIR/test.$OS.log"; # Log of results
  $self->{'continue_file'} = "$SCRIPT_DIR/script_restart_" . $self->arch;
  $self->setup_tests;
  return( $self );
}

sub process_args {
#
# Process the input arguments to the script
#
  my $self = shift;

  print "process_args:: Process the input arguments\n";
  while ( defined($ARGV[0]) ) {
    $_ = shift( @ARGV );
    if ( /^-help$/ ) {              # Help message
      $self->usage;
    # Clean the old directories out before running
    } elsif ( /^-noclean$/ ) {
      $self->do_clean( "no" );
    # If you want to continue a previous run of the script
    } elsif ( /^-r$|^-resume$/ ) {
      $self->do_continue( "yes" );
    # If you give the error-growth option use the base-line file from another machine
    } elsif ( /^-e$|^-errgro$/ ) {
      $_ = shift( @ARGV );
      if ( ! defined($_) ) { $self->usage; }
      # check for machine name with parenthesis around platform name:
      # And optionally, a :lab-name to specify the labname at remote site
      # for example babyblue(aix:ncar)
      if ( ! /(.+)\((.+)\)/ ) { $self->usage; }
      $self->{'remote_machine'}  = $1;
      $_ = $2;
      if ( ! /^([^:]+):*(.*)$/ ) { $self->usage; }
      $self->{'remote_platform'} = $1;
      $_ = $2;
      if ( /./  ) {
        $self->{'remote_lab'} = $_;
      } else {
        $self->{'remote_lab'} = $self->env("LAB");
      }
      # Make sure valid platform
      foreach my $os ( $self->OS_list ) {
        if ( $self->{'remote_platform'} eq $os ) { $_ = "yes"; }
      }
      if ( $_ ne "yes" ) { $self->usage; }
      # Make sure valid lab
      foreach my $os ( $self->{cam_lab_obj}->LAB_list ) {
        if ( $self->{'remote_lab'} eq $os ) { $_ = "yes"; }
      }
      if ( $_ ne "yes" ) { $self->usage; }
    # If you want to skip to a given test (or just a set of tests) (or pick a resolution)
    } elsif ( /^-s$|^-skip$/ ) {
      $_ = shift( @ARGV );
      if ( ! defined($_) ) { $self->usage; }
      if ( ! /([^:]+):([0-9]+)-*([0-9]*):*([x0-9.]*)L*([0-9]*)/ ) { $self->usage; }
      $self->CFG->setcfg( "DYNAMICS", $1 );  
      my $DYNAMICS = $self->CFG->cfg( "DYNAMICS", $1 );
      $self->{dynamics} = $1;
      my $test = $2; 
#
# Set resolution
#
      my $run_res = $self->{run_res};
      my $run_lev = $self->{run_lev};
      if ( $4 ne "" ) { 
         if ( $DYNAMICS ne "all" ) {
           $$run_res{$DYNAMICS} = $4; 
         } else {
           foreach my $dyn ( sort( keys(%$run_res) ) ) {
             $$run_res{$dyn} = $4; 
           }
         }
      }
      if ( $5 ne "" ) { 
         if ( $DYNAMICS ne "all" ) {
           $$run_lev{$DYNAMICS} = $5; 
         } else {
           foreach my $dyn ( sort( keys(%$run_res) ) ) {
             $$run_lev{$dyn} = $5; 
           }
         }
      }
      $self->{'firsttest'} = $test;
      $test = $test - 1;
      $self->{'lasttest'} = $3;
      if ( $self->{'lasttest'} eq "" ) { 
        $self->{'lasttest'} = undef; 
      } else {
        if ( $self->{'lasttest'} <= $test ) { $self->usage; }
      }
    # If you want to set the no fail option and continue even if something goes wrong
    } elsif ( /^-nofail$/ ) {
      $self->do_nofail( "yes" );
    # Compare to this model version
    } elsif ( /^-c$|^-compare$/ ) {
      $_ = shift( @ARGV );
      if ( ! defined($_) ) { $self->usage; }
      $self->do_compare( "yes" ); 
      $self->{'cam_cont_dir'} = $_;
      $self->{"CONT_ROOTDIR"} = $_;
    } elsif ( /^-l$|^-lab$/ ) {              # Set the Lab running at
      $_ = shift( @ARGV );
      if ( ! defined($_) ) { $self->usage; }
      $self->setenv( "LAB", $_ );     
    } elsif ( /^-fv2d$/ ) {              # Set the 2d parallel decomp for FV
      $self->{'twod_yz'} = "TRUE";
    } elsif ( /^-unique_id$/ ) {
      $_ = shift( @ARGV );
      if ( ! defined($_) ) { $self->usage; }
      $self->{'use_unique_id'} = "yes";
      my $unique_id = $_;
      $self->{"unique_id"} = $unique_id;
      # reset file paths and ensure unique names
      my $SCRIPT_DIR = $self->env( "SCRIPT_DIR" ); 
      my $CASE_DIR = $self->env( "CASE_DIR" ); 
      my $OS = $self->{'os'};
      # Log of results
      $self->{'testlog'} = "$SCRIPT_DIR/test.$OS.$unique_id.log";
      $self->{'continue_file'} = "$CASE_DIR/script_restart_" . $self->arch;
      $self->config_file("config.".$self->Platform.".$unique_id.cache.xml");
# BEGIN DEBUG
#      print "DEBUG:  ========================================\n";
#      print "DEBUG:  GOT OPTION -unique_id $_\n";
#      print "DEBUG:  use_unique_id = ".$self->{'use_unique_id'}."\n";
#      print "DEBUG:  unique_id = ".$self->{'unique_id'}."\n";
#      print "DEBUG:  testlog = ".$self->{'testlog'}."\n";
#      print "DEBUG:  continue_file = ".$self->{'continue_file'}."\n";
#      print "DEBUG:  config_file = ".$self->config_file."\n";
#      print "DEBUG:  EARLY EXIT FOR DEBUGGING\n";
#      print "DEBUG:  ========================================\n\n";
#      $self->usage;
# END DEBUG
    } else { 
      print "Invalid option $_\n\n";
      $self->usage;
    }
  }
  # If Env variable $CONT_ROOTDIR defined use it for comparision
  my $CONT_ROOTDIR = $self->{"CONT_ROOTDIR"};
  if ( ! $self->do_compare && defined($CONT_ROOTDIR) ) {
    $self->do_compare( "yes" ); 
    $self->{'cam_cont_dir'} = $CONT_ROOTDIR;
  }
  $self->do_build( "yes" ); 
  $self->do_build_namelist( "yes" ); 
  $self->do_log( "yes" ); 

  print "Build and run test simulations\n";
}

sub usage {
#
# Usage statement if command arguments are done correctly
#
  my $self = shift;

  my @dyn = $self->dynamics;
  my @os_list = $self->OS_list;
  my @LAB_list = $self->{cam_lab_obj}->LAB_list;
  my $LAB = $self->env( "LAB" );
  print <<EOF;
Usage: perl $0 [options]

Options are:

	-lab     = Set the lab you are running at 
                 (of @LAB_list) [$LAB]
                 (Also set by setting the env variable LAB)
	-help    = Help (this message, also lists tests performed)
	-noclean = Don't clean the old directories out.
	-nofail  = Continue even if errors are found
	-resume  = Continue a previous run of the script at the point it left 
	           off at.  If used with -unique_id, see caveats in the
		   description of that option.  
	-errgro mach(plat:lab) = List the remote machine and platform to use as the 
                 baseline for error-growth tests (ie. -e "babyblue.ucar.edu(aix:ncar)" )
                 (list of valid platforms are: @os_list)
                 (tests must be in the default location on the remote machine)
                 (list of valid labs are: @LAB_list [leave off :lab to use $LAB])
	-fv2d    = Turn on the 2D parallel decomposition for FV dynamics.
	-skip dy:#:res = Skip to given dynamics and test # (or range of numbers)
                 (example -skip sld:9 start with sld dynamics test no. 9)
                 (or      -skip fv:2-4 start with fv dynamics test no. 2 and do up to test 4)
                 (or      -skip all:2 start at tests 2 for all dynamics)
                 (or      -skip eul:2-4:64x128L26 do tests 2 through 4 for eul at 64x128 (T42) with 26 levs)
	-compare cam-root-dir = Compare to given version of the model in this directory
                 This is the root directory of the cam1 directory tree 
                        (Example -compare /home/erik/cam1)
                        (Also set by setting the env variable CONT_ROOTDIR)
	-unique_id job_id    Use this if you want to run more than one 
	         instance of test-model concurrently.  If different values are 
		 supplied for each job_id, then the instances should be able 
		 to run independently.  If -resume is used to restart a 
		 previous test, then job_id must be set to the same value used 
		 in the previous test.  

Batch operation: To submit to batch you need to create and/or edit the
      submission information at the top of a script. This directory contains
      simple batch submission scripts for each lab that can be used for this
      purpose. Most likely the que name and possibly the number of nodes might 
      need to be changed. You need to edit either the PBS, NQS, or Loadleveler 
      sections, depending on which queing system your machine supports.

      To specify non-standard options for batch operation edit the execution
      statement in the batch script for $0.

      To submit script ${LAB}_batch.csh to the batch queues, follow the 
      instructions in the script.  

EOF
  $self->list_tests;   # List the tests out
  die "Terminating";
}

sub setup_tests {
#
# Setup the tests that will be run
# If you change the tests here you also may need to change test-model.pl
#
  my $self = shift;

  my %run_res = (eul=>"32x64", sld=>"32x64", fv=>"4x5");# Horizontal resolution
  $self->{run_res} = \%run_res;
  my %run_lev = (eul=>26,    sld=>26,    fv=>26);   # Number of vertical levels
  $self->{run_lev} = \%run_lev;
#
# Tests to run and run-type for each, also give the number of steps to run for each
# (The tests will run in alphabetical order, so include a letter to order it)
#
# Save/restart tests should be timed to maximize the number of restart files
# that need to be saved. Hence, odd times should be used, so that both abs/ems
# restart datasets are saved as well as history restart datasets.
# (If you change the keys you need to change the script later as well)
#
# Description of tests should start with a number and then have a "_".
# They also should be in order.
#
  my %debug_config= ( '01_initial_debug_run_SPMD'     =>{ SPMD=>"TRUE",  SMP=>"TRUE" },
                      '02_initial_debug_run_SPMDnoSMD'=>{ SPMD=>"TRUE",  SMP=>"FALSE" },
                      '03_initial_debug_run_nonSPMD'  =>{ SPMD=>"FALSE", SMP=>"TRUE" } );
  $self->{debug_config} = \%debug_config;
# Run-tests for full-physics
  my %runtypes = ( 
      '04_initial'=>"initial", 
      '05_restart'=>"restart",
      '06_initial_compare_to_restart'=>"initial" );
  $self->{runtypes} = \%runtypes;
  my %nestep  = ( 
      eul=>{'04_initial'=>5, '05_restart'=>10, '06_initial_compare_to_restart'=>30},
      sld=>{'04_initial'=>5, '05_restart'=>10, '06_initial_compare_to_restart'=>10},
       fv=>{'04_initial'=>5, '05_restart'=>10, '06_initial_compare_to_restart'=>10} );
  $self->{nestep} = \%nestep;
# Run important physics test
  my %physics_tests = (
        "07_SOM"=>{OCEANMODEL=>"som"},
       );
  $self->{physics_tests} = \%physics_tests;
# Run non-perturbation test with control code
  my %control_tests = (
        "08_control_nonpert"=>{OCEANMODEL=>"dom"},
        "09_control_SOM"=>{OCEANMODEL=>"som"},
       );
  $self->{control_tests} = \%control_tests;
# Error-growth run-tests
  $self->{roundoff}     = 1e-14;
  my $roundoff = $self->{roundoff};
  my %errgrotests  = ( 
      '10_error_growth_adiabatic'=>{adiabatic=>"true",pertlim=>0},
      '11_error_growth_adiabatic_pert'=>{adiabatic=>"true",pertlim=>$roundoff},
      '12_error_growth_full_physics'=>{adiabatic=>"false",pertlim=>0},
      '13_error_growth_full_physics_pert'=>{adiabatic=>"false",pertlim=>$roundoff}
  );
  $self->{'errgrotests'} = \%errgrotests;
# Control-run tests
  my %controltests = ( 
       '14_control_adiabatic'=>{adiabatic=>"true", pertlim=>0},
       '15_control_adiabatic_pert'=>{adiabatic=>"true", pertlim=>$roundoff},
       '16_control_full_physics'=>{adiabatic=>"false", pertlim=>0},
       '17_control_full_physics_pert'=>{adiabatic=>"false", pertlim=>$roundoff} );
  $self->{controltests} = \%controltests;
# Special tests (only ran under specific dynamics)
  my %specialtests = ( 
       '18_dataicemodel'=>{DYNAMICS=>"eul", LANDMODEL=>"CLM2", PHYSICS=>"cam1", 
                              OCEANMODEL=>"dom", ICEMODEL=>"ccmice", RESOLUTION=>undef, 
                              PLEV=>undef},
     );
  $self->{specialtests} = \%specialtests;
}

sub dynamics {
#
# List dynamics that will be gone through
#
  my $self = shift;

  my @dyn;
  if ( $self->{'dynamics'} eq "all" ) {
    my $run_res = $self->{run_res};
    @dyn = sort( keys( %$run_res ) );
  } else {
    @dyn = ($self->{'dynamics'});
  }
  return( @dyn );
}

sub tests {
#
# List tests that will be gone through for a given class of tests
#
  my $self = shift;
  $_ = shift;

  my $nm = ref($self) . "::tests";
  my @list = ();
  SWITCH: {
    (/debug/) && do {
         my $debug_config= $self->{debug_config};
         @list = sort( keys( %$debug_config) );
       last SWITCH;
    };
    (/restart/) && do {
         my $runtypes = $self->{runtypes};
         @list = sort( keys( %$runtypes) );
       last SWITCH;
    };
    (/physics/) && do {
         my $physics = $self->{physics_tests};
         @list = sort( keys( %$physics) );
       last SWITCH;
    };
    (/control/) && do {
         if ( $self->do_compare ) {
           my $control = $self->{control_tests}; 
           @list = sort( keys( %$control) );
         } else {
           @list = ();
         }
       last SWITCH;
    };
    (/^errgro/) && do {
         my $errgrotests = $self->{'errgrotests'};
         foreach my $test ( sort( keys( %$errgrotests ) ) ) {
           if ( defined($$errgrotests{$test}{'adiabatic'}) ) {
             push( @list, $test );
           }
         }
       last SWITCH;
    };
    (/conterrgro/) && do {
         if ( $self->do_compare ) {
           my $controltests = $self->{controltests};
           @list = sort( keys( %$controltests ) );
         } else {
           @list = ();
         }
       last SWITCH;
    };
    (/special/) && do {
         my $specialtests = $self->{specialtests};
         @list = sort( keys( %$specialtests ) );
       last SWITCH;
    };
#---------------------------------------------------------------------------------------
# default
#---------------------------------------------------------------------------------------
    die "ERROR::($nm) bad option sent to tests\n";
  }
#
# Loop over tests and make sure within bounds specified on command-line
#
  my @newlist = ();
  my $test;
  while ( defined($test = shift(@list)) ) {
     if ( $test !~ /(^\d+)_/ ) {
       die "ERROR::($nm) bad test name\n";
     }
     my $testno = $1;
     if ( ! defined($self->{'firsttest'}) || ! defined($self->{'lasttest'}) 
     || (($testno >= $self->{'firsttest'}) && ($testno <= $self->{'lasttest'})) ) {
       push( @newlist, $test );
     }
  }
  return(@newlist);
}

sub config {
#
# Return given configuration setting for a given test
#
  my $self = shift;
  my $type = shift;
  my $desc = shift;

  my $nm = ref($self) . "_config";
  my $value;
  if ( ! defined($type) ) {
    die "ERROR::($nm):: type not sent to config: $type\n";
  }
  my $dyn = $self->CFG->cfg("DYNAMICS");
  if (    $type eq "TWOD_YZ" ) {
    if ( $dyn eq "fv" ) {
      $value = $self->{twod_yz};
    } else {
      $value = "FALSE";
    }
  } elsif ( ! defined($desc) ) {
    if (    $type eq "PLEV" ) {
      my $run_lev = $self->{run_lev};
      $value = $$run_lev{$dyn};
    } elsif ( $type eq "RESOLUTION" ) {
      my $run_res = $self->{run_res};
      $value = $$run_res{$dyn};
    } elsif ( $type eq "roundoff" ) {
      $value = $self->{roundoff};
    } else {
      die "ERROR::($nm):: Invalid type sent to config: $type\n";
    }
  } else {
    if (    $type eq "RUNTYPE" ) {
      my $runtypes = $self->{runtypes};
      $value = $$runtypes{$desc};
    } elsif ( $type eq "SPMD" ) {
      my $debug_config= $self->{debug_config};
      $value = $$debug_config{$desc}{SPMD};
    } elsif ( $type eq "SMP" ) {
      my $debug_config= $self->{debug_config};
      $value = $$debug_config{$desc}{SMP};
    } elsif ( $type eq "nestep" ) {
      my $nestep = $self->{nestep};
      $value = $$nestep{$dyn}{$desc};
    } elsif ( $type eq "adiabatic" ) {
      my $errgrotests = $self->{'errgrotests'};
      $value = $$errgrotests{$desc}{'adiabatic'};
      if ( ! defined($value) ) {
        my $controltests = $self->{controltests};
        $value = $$controltests{$desc}{'adiabatic'};
        if ( ! defined($value) ) {
          die "ERROR::($nm):: $desc test does not define adiabatic value\n";
        }
      }
    } elsif ( $type eq "pertlim" ) {
      my $errgrotests = $self->{'errgrotests'};
      $value = $$errgrotests{$desc}{'pertlim'};
      if ( ! defined($value) ) {
        my $controltests = $self->{controltests};
        $value = $$controltests{$desc}{'pertlim'};
        if ( ! defined($value) ) {
          die "ERROR::($nm):: $desc test does not define pertlim value\n";
        }
      }
    } elsif ( $type eq "PLEV" ) {
      my $specialtests = $self->{specialtests};
      $value = $$specialtests{$desc}{PLEV};
    } elsif ( $type eq "DYNAMICS" ) {
      my $specialtests = $self->{specialtests};
      $value = $$specialtests{$desc}{DYNAMICS};
    } elsif ( $type eq "LANDMODEL" ) {
      my $specialtests = $self->{specialtests};
      $value = $$specialtests{$desc}{LANDMODEL};
    } elsif ( $type eq "ICEMODEL" ) {
      my $specialtests = $self->{specialtests};
      $value = $$specialtests{$desc}{ICEMODEL};
    } elsif ( $type eq "OCEANMODEL" ) {
      if ( $desc =~ /_control/ ) {
        my $control_tests = $self->{control_tests};
        $value = $$control_tests{$desc}{OCEANMODEL};
      } elsif ( $desc =~ /[0-9]+_SOM/ ) {
        my $physics_tests = $self->{physics_tests};
        $value = $$physics_tests{$desc}{OCEANMODEL};
      } else {
        my $specialtests = $self->{specialtests};
        $value = $$specialtests{$desc}{OCEANMODEL};
      }
    } elsif ( $type eq "RESOLUTION" ) {
      my $specialtests = $self->{specialtests};
      $value = $$specialtests{$desc}{RESOLUTION};
    } elsif ( $type eq "PHYSICS" ) {
      my $specialtests = $self->{specialtests};
      $value = $$specialtests{$desc}{PHYSICS};
    } else {
      die "ERROR::($nm):: Invalid type sent to config: $type, description=$desc\n";
    }
  }
  return( $value );
}

sub list_tests {
#
# List the tests that are performed.
#
  my $self = shift;

  my @list = $self->tests("debug");
  push @list, $self->tests("restart");
  push @list, $self->tests("physics");
  push @list, $self->tests("control");
  push @list, $self->tests("errgro");
  push @list, $self->tests("conterrgro");
  print <<EOF;
  List of tests that are performed:

  for each dynamics at these resolutions:

EOF
  my $run_res = $self->{run_res};
  my $run_lev = $self->{run_lev};
  foreach my $dyn ( $self->dynamics ) {
    print "$dyn: \tHorizontal Resolution: $$run_res{$dyn} " . 
          "# of vertical levels: $$run_lev{$dyn}\n";
  }
  foreach my $i ( @list ) {
    print "$i\n";
  }
  my @splist = $self->tests("special");
  my $specialtests = $self->{specialtests};
  foreach my $i ( @splist ) {
    print "$i only " . $$specialtests{$i}{DYNAMICS} . " dynamics\n";
  }
  die "\n\nTerminating";
}

sub do_remote_errg {
#
# Test if using a remote error-growth file as the baseline
#
  my $self = shift;
  my $mach = $self->{'remote_machine'};
  my $plat = $self->{'remote_platform'};
  if ( defined($mach) && defined($plat) ) {
    return( 1 );
  } else {
    return( 0 );
  }
}

sub remote_baseline_errg {
#
# Return Name of the remote error-growth filename to compare to
# after copying it (using scp) from the remote to local location.
# scp is assumed to be in your path, or else a error will occur.
#
  my $self = shift;
  my $file = shift;

  my $mach = $self->{'remote_machine'};

  my $local_filename;
  if ( ! $self->do_remote_errg ) {
    $local_filename = $file;
  } else {
    # Seperate directory name from filename
    my @path = split( "/", $file );
    my $filename = $path[$#path];
    # Figure out remote pathname of file to get
    # Swap platform name of current case with platform name of remote
    my $rem = $self->duplicate;
    $rem->setenv( "LAB", $self->{'remote_lab'} );
    $rem->OS( $self->{'remote_platform'} );
    my $remote_case_dir = $rem->{cam_lab_obj}->default( "CASE_DIR", $rem );
    my $remote_plat = $rem->Platform;
    my $remote_case = $rem->env( "CASE" );
    my $local_plat = $self->Platform;
    $remote_case =~ s/$local_plat/$remote_plat/;
    my $remote_dir = "$remote_case_dir/$remote_case";
    my $remote_filename = "$mach:$remote_dir/$filename";
    $local_filename = "$file.$mach.".$rem->OS;
    # Copy file here
    $self->exec( "scp $remote_filename $local_filename", 
                 "Error:: copying remote baseline file" );
  }
  return( $local_filename );
}

sub setup_directories {
#
# Extend the set_directories method to handle comparing to previous
# versions.
#
  my $self = shift;
  my $comp = shift;

  if ( defined($comp) && ($comp eq "compare") ) {
    print "Setup directories for comparision: ";
    if ( defined($self->{'cam_cont_dir'}) ) {
      #
      # Change CAMROOT to comparision directory, change back when continuation_mark called
      #
      my $CAMROOT = $self->CFG->cfg( "CAMROOT" );
      my $MODEL_CFGDIR = $self->CFG->cfg( "MODEL_CFGDIR" );
      $self->{'rootdir'} = $CAMROOT;
      if ( ! defined($MODEL_CFGDIR) ) {
        $MODEL_CFGDIR = "$CAMROOT/models/atm/cam/bld";
      }
      $self->{'cfgdir'}  = $MODEL_CFGDIR;
      $CAMROOT = $self->{'cam_cont_dir'};
      $self->CFG->setcfg( "CAMROOT", $CAMROOT );
      $MODEL_CFGDIR = "$CAMROOT/models/atm/cam/bld";
      $self->CFG->setcfg( "MODEL_CFGDIR", $MODEL_CFGDIR );
      print "compare to $CAMROOT\n";
    } else {
      print "\n";
    }
    $self->CFG->unsetcfg( "MODEL_BLDDIR" );
    $self->CFG->unsetcfg( "MODEL_EXEDIR" );
    $self->SUPER::setup_directories;
  } elsif ( defined($comp) && ($comp eq "support") ) {
    print "Setup support directories: \n";
    #
    # Build the Case directory (where to build code and run cases)
    #
    my $OS = $self->OS;
    my $CASE_DIR = $self->env( "CASE_DIR" );
    if ( ! defined($CASE_DIR) ) {
      $CASE_DIR = $self->{cam_lab_obj}->default( "CASE_DIR", $self );
    }
    if ( ! -d $CASE_DIR ) {
       mkdir( $CASE_DIR, 0755 ) || die "Can not mkdir $CASE_DIR";
    }
    print "Case directory: $CASE_DIR\n";
    $self->setenv( "CASE_DIR", $CASE_DIR );
    my $BUILD_DIR = $self->env( "BUILD_DIR" );
    if ( ! defined($BUILD_DIR) ) {
      $BUILD_DIR = $self->{cam_lab_obj}->default( "BUILD_DIR", $self );
    }
    if ( ! -d $BUILD_DIR ) {
      mkdir( $BUILD_DIR, 0755 ) || die "Can not mkdir $BUILD_DIR";
    }
    $self->setenv( "BUILD_DIR", $BUILD_DIR );
    my $dir = "$BUILD_DIR/cprncobj";
    if ( ! -d $dir ) {
      mkdir( $dir, 0755 ) || die "Can not mkdir $dir";
    }
    system( "/bin/rm $dir/Makefile" );
    my $MODEL_CPRDIR = $self->env( "MODEL_CPRDIR" );
    symlink( "$MODEL_CPRDIR/Makefile", "$dir/Makefile" );
  } else {
    $self->SUPER::setup_directories;
  }
}

sub do_nofail {
#
#  Test if the no fail option is set
#
  my $self = shift;
  my $value = shift;

  return( $self->do('nofail', $value) );
}

sub do_compare {
#
#  Test if should compare to a control source directory or not
#
  my $self = shift;
  my $value = shift;

  return( $self->do('compare', $value) );
}

sub do_continue {
#
#  Test if should continue a previous run of the script or not
#
  my $self = shift;
  my $value = shift;

  return( $self->do('continue', $value) );
}

sub completed {
#
#  Test if this part has already been completed or not
#
  my $self = shift;
  my $value = shift;

  #
  # If just begining mark completed as false so that 
  # everything will be done until it gets to the point of
  # running tests.
  #
  if ( $self->do('begin') ) {
    return( 0 );
  # Now return completed based on the value of COMPLETED 
  } else {
    return( $self->do('completed', $value) );
  }
}

sub numerically {
#
# Use with "sort" to sort values in numerical order
#
  $a <=> $b;
}


sub build_namelist {
#
# Build namelist only if you aren't continuing.
#
  my $self = shift;
  my $reset = shift;

  if ( $self->do_continue && $self->completed ) {
    print "Continue option set, skip build_namelist\n";
    return;
  }
  $self->SUPER::build_namelist( $reset );
}

sub exec {
#
# Extend the exec method so that if continuation option is set and this
# section already done, it doesn't actually execute.
#
  my $self = shift;
  my $command = shift;
  my $die_msg = shift;
  my $echo = shift;

  if ( $self->do_continue && $self->completed ) {
    print "Continue option set, skip running: $command\n";
    return;
  }
  $self->SUPER::exec( $command, $die_msg, $echo );
}

sub run {
#
# Extend the run method so that if continuation option is set and this
# section already done, it doesn't actually do anything.
#
  my $self = shift;
  my $desc = shift;
  my $logfile = shift;


  if ( $self->do_continue && $self->completed ) {
    print "Continue option set, skip running model: $desc\n";
    if ( defined($logfile) ) { $self->{'logfile'} = $logfile; }
    return;
  }
  $self->SUPER::run( $desc, $logfile );
}


sub compare_files {
#
# Compare the history files using "cprnc"
#
  my $self = shift;
  my $file1 = shift;
  my $file2 = shift;
  my $log = shift;
  my $die_msg = shift;

  print "Compare the two history tapes $file1 and $file2 output to $log:\n";
  if ( $self->do_continue && $self->completed ) {
    print "Continue option set, skip this file comparision\n";
    return;
  }
  my $CASE_DIR = $self->env( "CASE_DIR" );
  $self->exec( "$CASE_DIR/cprnc $file1 $file2 > $log", "Error:: running cprnc" );
  open( LOG, "<$log" ) || $self->die( "Error, opening cprout file: $log" );
  my $different = 0;
  while( defined($_ = <LOG>) ) {
    if ( /RMS ([a-zA-Z0-9_-]+) (.+)/ ) {
      $CAM_test::cprout{$log}{$1} = $2;
      if ( $2 !~ /0.[0]+[eE]\+[0]+/ ) { 
        $different = 1; 
      }
    }
  }
  close( LOG );
  if ( defined($die_msg) ) {
    if ( $different ) {
      $self->die( "Error: $die_msg, see $log" );
    } else {
      print "RMS difference between all fields on the history files are identical!" .
            " Comparision successful, continuing...\n";
    }
  }
  return( $different );
}

sub report_compare {
#
# Report a summary of the comparisions
#
  my $self = shift;
  my @cprlogs = @_;

  my $log; my $var; my %diff;
  #
  # Loop through comparision logs and report which fields bit-for-bit or different
  #
  foreach $log ( @cprlogs ) {
    my $ref = $CAM_test::cprout{$log};
    if ( defined($ref) ) {
      print "For $log:\n";
      my %var_diffs = %$ref;
      my @var_list = sort( keys(%var_diffs) );
      my @bit4bit_list;
      my @diff_list;
      my $max_diff = 0.0;
      foreach $var ( sort( keys(%var_diffs) ) ) {
        if ( $var_diffs{$var} > $max_diff ) { $max_diff = $var_diffs{$var}; }
        if ( $var_diffs{$var} =~ /0.[0]+[eE]\+[0]+/ ) {
          push @bit4bit_list, $var;
        } else {
          push @diff_list, $var;
        }
      }
      print "Fields on files: @var_list\n";
      $diff{$log} = "NOT";
      if ( $#bit4bit_list == $#var_list ) { 
        print "All fields are bit for bit\n"; 
        $diff{$log} = "b4b";
      } elsif ( $#diff_list == $#var_list ) { 
        print "No fields are bit for bit\n"; 
        print "Max difference: $max_diff\n";
      } else {
        print $#bit4bit_list+1 . " fields are bit for bit\n";
        print $#diff_list+1 . " fields are different: @diff_list\n";
        print "Max difference: $max_diff\n";
      }
    }
  }
  my $string;
  #
  # Report on success of testing
  #
  if ( $self->{'errors'} == 0 ) {
    $string = "\n\nModel testing was successful!\n";
  } else {
    $string = "\n\nERROR!!!:: ".$self->{'errors'}." errors found during execution!\n";
  }
  my $testlog = $self->{'testlog'};
  open( TESTLOG, ">>$testlog" ) || die "Can not open $testlog";
  print TESTLOG $string;
  #
  # Now Loop through logs and report if bit-for-bit with control library
  #
  my $b4b = "yes";
  foreach $log ( @cprlogs ) {
    my $ref = $CAM_test::cprout{$log};
    if ( defined($ref) && ($log =~ /con.+control.+/) && ($log !~ /_pert/) ) {
      if ( $diff{$log} eq "NOT" ) { $b4b = "NO"; }
    }
  }
  if ( $self->do_compare ) {
    if ( $b4b eq "yes" ) {
      $string = "Code is bit-for-bit with control library\n";
    } else {
      $string = "\n\nWARNING!!!:: Code is NOT bit-for-bit with control library\n" .
            "             Verify that this is ok or validate that at least within roundoff.\n\n";
    }
    print TESTLOG $string;
  }
  #
  # Loop through error-growth logs and report if error-growth is high
  # (Just checking T)
  #
  foreach $log ( @cprlogs ) {
    my $ref = $CAM_test::cprout{$log};
    if ( defined($ref) && ($log =~ /errg.+error_growth.+/) ) {
      my %var_diffs = %$ref;
      if ( $var_diffs{T} > 1.e-5 ) { 
        $string = "\n\nWARNING!!!:: Error growth seems to be high!\n" .
              "             Verify that this is ok.\n\n";
        print TESTLOG $string;
      }
    }
  }
  #
  # Report on performance comparision
  #
  my $timing = $self->{'timing'};
  my %comp = ( '06_initial_compare_to_restart'=>"08_control_nonpert",
               '10_error_growth_adiabatic'=>"14_control_adiabatic",
               '12_error_growth_full_physics'=>"16_control_full_physics" );
  foreach my $dyn ( $self->dynamics ) {
     foreach my $test ( keys(%comp) ) {
       if ( defined($$timing{$dyn}{$test}) &&
            defined($$timing{$dyn}{$comp{$test}}) ) {
         my $rate = $$timing{$dyn}{$test}->compare_perf( $$timing{$dyn}{$comp{$test}} );
         if ( $rate < -15.0 ) {
           my $string = "\n\nWARNING!!!:: Model is significantly slower than control by ", 
                        -$rate, 
                        "\n               comparing case: $test to $comp{$test}!\n";
           print TESTLOG $string;
         }
       }
     }
  }
  close( TESTLOG );
}

sub lower_PEs {
#
# Change the number of PE's for testing
#
  my $self = shift;

  print "Lower the number of PE\'s used to ensure that different configurations match results:\n";
  my $SHMEM_CPUS = $self->env( "SHMEM_CPUS" );
  if ( defined($SHMEM_CPUS) ) {
    $CAM_run::SAVE_SHMEM_CPUS = "$SHMEM_CPUS";   # Save number it was set to
    print "Current number of shared memory CPUS: " . $SHMEM_CPUS. "\n";
    $SHMEM_CPUS = $SHMEM_CPUS / 2;
    if ( $SHMEM_CPUS < 1 ) { 
      $SHMEM_CPUS = 1; 
      print "Leave number at:                      " . $SHMEM_CPUS. "\n";
    } else {
      print "Change number to:                     " . $SHMEM_CPUS. "\n";
    }
  } else {
    $CAM_run::SAVE_SHMEM_CPUS = 0;
    $SHMEM_CPUS = 1;
    print "Change number of shared memory CPUS:  " . $SHMEM_CPUS. "\n";
  }
  $self->setenv( "SHMEM_CPUS", $SHMEM_CPUS );
  my $SPMD_NODES = $self->env( "SPMD_NODES" );
  my $SPMD = $self->CFG->cfg( "SPMD" );
  if ( ($SPMD eq "TRUE") && defined($SPMD_NODES) ) {
    $CAM_run::SAVE_SPMD_NODES = "$SPMD_NODES";   # Save number it was set to
    print "Current number of SPMD nodes: " . $SPMD_NODES . "\n";
    $SPMD_NODES = $SPMD_NODES / 2;
    if ( $SPMD_NODES < 2 ) {
      $SPMD_NODES = 2;
      if ( ! defined($SHMEM_CPUS) ) {
        $CAM_run::SAVE_SHMEM_CPUS = 0;
      }
      $SHMEM_CPUS = 1;
      print "Keep SPMD nodes at $SPMD_NODES change shared-memory CPU\'s to $SHMEM_CPUS\n";
    } else {
      print "Change number to            : " . $SPMD_NODES . "\n";
    }
    $self->setenv( "SPMD_RUNCMND", undef ); # un-define the SPMD run command so that node change will take effect
                           # this happens when run_time_env method is invoked
    $self->setenv( "SPMD_NODES", $SPMD_NODES );
  }
}

sub start_log {
#
# Start the log-file
#
  my $self = shift;

  my $testlog = $self->{'testlog'};
  if ( -f "$testlog" ) {
    `/bin/mv $testlog $testlog.last`;
  }
  open( TESTLOG, ">$testlog" ) || die "Can not open $testlog";
  my $hostname = `/bin/hostname`; chomp( $hostname );
  my $date = `date`;
  my $OS = $self->OS;
  print TESTLOG "Log of test results for $hostname a $OS on $date\n\n";
  close( TESTLOG );
}

sub read_mark {
#
# Read the restart filemark
#
  my $self = shift;

  my $SCRIPT_DIR = $self->env( "SCRIPT_DIR" );
  my $filename = $self->{continue_file};
  open( MARK, "<$filename" ) || die "Could not open $filename";
  $_ = <MARK>;
  my $read_mark;
  if ( ! /(^[A-Za-z: ]+) ([0-9]+_)/ ) {
    die "ERROR:: Continuation mark on $filename is not in correct format\n";
  }
  my $dyn  = $1;
  my $desc = $2;
  if ( $dyn =~ /all/ ) {
    $read_mark = $self->mark( $desc );
  } else {
    $read_mark = "$dyn $desc";
  }
  close( MARK );
  return( $read_mark );
}

sub write_mark {
#
# Write the restart file-mark out
#
  my $self = shift;
  my $desc = shift;

  my $SCRIPT_DIR = $self->env( "SCRIPT_DIR" );
  my $filename = $self->{continue_file};
  my $mark = $self->mark( $desc );
  open( MARK, ">$filename" ) || die "Could not open $filename";
  print MARK "$mark\n";
  close( MARK );
}

sub mark {
#
# Return mark to designate how far the testing has proceeded
# description assumed to be a integer number followed by an optional description.
#
  my $self = shift;
  my $desc = shift;

  # Check that desc set
  if ( ! defined($desc) ) { die "Error: desc not sent to mark method"; }
  # Make sure begining number printed out as two characters
  my $number = $desc;
  my $rest;
  if ( $desc =~ /(^[0-9]+)(.*)$/ ) {
    $number = sprintf "%2.2d", $1;
    $rest   = $2;
  }
  my $DYNAMICS = $self->CFG->cfg( "DYNAMICS" );
  return( "Dynamics: $DYNAMICS $number$rest" );
}

sub continuation_mark {
#
# Mark and/or check how far the script has completed.
#
  my $self = shift;
  my $desc = shift;

  if ( defined($self->{'rootdir'}) ) {
    $self->CFG->setcfg( "CAMROOT", $self->{'rootdir'} );
    $self->CFG->setcfg( "MODEL_CFGDIR", $self->{'cfgdir'} );
    $self->{'rootdir'} = undef;
    $self->{'cfgdir'} = undef;
  }
  my $SCRIPT_DIR = $self->env( "SCRIPT_DIR" );
  my $DYNAMICS = $self->CFG->cfg( "DYNAMICS" );
  my $filename = $self->{continue_file};
  my $mark = $self->mark( $desc );
# TBH:  Moved this to the end after discussion with Erik to fix bug that 
# TBH:  caused "ran" to be printed even when a test failed.  
#  $self->{'status'} = "ran";
  # If script is completed
  if ( $desc eq "done" ) {
    `/bin/rm $filename`;
    my $testlog = $self->{'testlog'};
    print "Results of the tests are in $testlog\n";
    open( TESTLOG, "<$testlog" ) || die "Can not open $testlog";
    while( defined($_ = <TESTLOG>) ) {
      print $_;
    }
    close( TESTLOG );
  # If script is just starting
  } elsif ( ($desc =~ /^0_/) && ($DYNAMICS eq "none") ) {
    $self->do('begin',"no");
    $self->start_log;
  # If you are continuing and haven't found the part that isn't done yet
  } elsif ( $self->do_continue && $self->completed ) {
    my $read = $self->read_mark;
    # Set completed to no if reached the section where it stopped
    if ( $mark =~ /^$read/ ) {
      $self->completed( "no" );
    }
    $self->{'status'} = "skipped";
  } else {
    $self->write_mark( $desc );
    #
    # Print test results to log file
    #
    if ( $desc !~ /^0_/) {
      my $RES = $self->CFG->cfg( "RESOLUTION" );
      my $LEV = $self->CFG->cfg( "PLEV" );
      my $testlog = $self->{'testlog'};
      open( TESTLOG, ">>$testlog" ) || die "Could not access testlog file: " . 
                                            $testlog;
      print TESTLOG  "$mark: " . "${RES}L${LEV} " . $self->{'status'} . "\n";
      close( TESTLOG );
    }
    #
    # If just doing part of the tests stop if reached the last test
    #
    my $lasttest = $self->{'lasttest'};
    if ( defined($lasttest) ) {
      my $last = $self->mark( $lasttest );
      $self->{'status'} = "ran";
      # If at last test mark completed to yes and skip the rest of the tests
      if ( $mark =~ /^$last/ ) {
        $self->completed( "yes" );
      }
    }
  }
  # TBH:  reset status to successful for next test
  # TBH:  assumption = continuation_mark is called after each test finishes 
  # TBH:  for either "success" or "failure" cases.  
  $self->{'status'} = "ran";
}

sub reset_PEs {
#
# Reset the number of PE's changed by "lower_PEs" to the original values
#
  my $self = shift;

  print "Reset the PE\'s\n";
  if ( defined($self->env("SPMD_RUNCMND")) ) {
    $self->setenv( "SPMD_RUNCMND", undef );
  }
  if ( defined($self->env("SHMEM_CPUS")) ) {
    $self->setenv( "SHMEM_CPUS", $CAM_run::SAVE_SHMEM_CPUS );
    if ( $CAM_run::SAVE_SHMEM_CPUS != 0 ) {
      print "Shared memory CPU\'s: ".$self->env( "SHMEM_CPUS" ). "\n";
    } else {
      print "Reset shared memory CPU\'s to default on machine:\n";
    }
  }
  if ( defined($self->env("SPMD_NODES")) ) {
    $self->setenv( "SPMD_NODES", $CAM_run::SAVE_SPMD_NODES );
    print "SPMD nodes: ".$self->env( "SPMD_NODES" ). "\n";
  }

}

sub die {
#
# Terminate gracifully if found an error
#
  my $self = shift;
  my $desc = shift;

  if ( $self->do_nofail ) {
    my $testlog = $self->{'testlog'};
    open( TESTLOG, ">>$testlog" ) || die "Could not access testlog file: " . 
                                          "$testlog";
    print TESTLOG "TEST FAILED -- $desc\n";
    close( TESTLOG );
    $self->{'status'} = "FAIL";
    $self->{'errors'} += 1;
  } else {
    die "$desc";
  }
}

sub version {
#
# Get the version from the ChangeLog
#
  my $self = shift;

  my $nm = ref($self) . "::version";
  my $ROOT_DIR =  $self->CFG->cfg( "CAMROOT" );
  my $file = "$ROOT_DIR/models/atm/cam/doc/". $self->{CHANGELOG};
  open( FILE, "<$file" ) || die "ERROR:($nm) Opening $file\n";
  $_ = <FILE>; 
  if ( ! defined($_) ) {
    die "ERROR::($nm) Reading file: $file\n";
  }
  if ( ! defined($_) ) {
    die "ERROR::($nm) Reading file: $file\n";
  }
  $_ = <FILE>;
  my $version = <FILE>;
  if ( ! defined($version) ) {
    die "ERROR::($nm) Reading file: $file\n";
  }
  chomp( $version );
  if ( $version !~ /cam|ccm[0-9]+_[0-9]+/ ) {
    die "ERROR::($nm) Reading version in file: $file\n";
  }
  return( $version );
}

sub save_timing {
#
# Save the timing information for this case
#
  my $self = shift;
  my $desc = shift;
  my $logfile = shift;

  my $DYNAMICS     = $self->CFG->cfg( "DYNAMICS" );
  my $MODEL_EXEDIR = $self->CFG->cfg( "MODEL_EXEDIR" );
  my $CASE = "$DYNAMICS:$desc";
  my $time = $self->{timing};
  if ( ! -f "$logfile" ) {
    return;
  }
  my %opt = ( dir=>"$MODEL_EXEDIR", log=>"$logfile", case=>"$CASE", nofail=>$self->do_nofail );
  $$time{$DYNAMICS}{$desc} = cam_timing->new( \%opt );
  my $return = $$time{$DYNAMICS}{$desc}->report_perf;
  if ( $self->do_nofail && ($return != 0) ) {
    $self->die( "ERROR!!!:: Error reporting on performance!\n" );
  }
}


1   # To make use or require happy
