#!/usr/bin/env perl
#=======================================================================
#
#  This is a CAM stand-alone atm model test script
#
# Usage:
#
# perl test-model.pl options
#
#=======================================================================

use 5.004;   # Use at least this version of perl
use Cwd;     # Use current working directory module
#
# Use the shortened syntax for the following ENV variables
#
use Env qw(ICEMODEL LANDMODEL OCEANMODEL PHYSICS SCRIPT_DIR SPMD SMP SPMD_NODES F_OPTIMIZATION_OVERRIDE USER_FC);
use strict;
#use diagnostics;

# SCRIPT_DIR location of the main script
if ( defined($SCRIPT_DIR) ) {
  chdir( $SCRIPT_DIR );
} else {
  $SCRIPT_DIR = cwd( );
}

use lib ".", "../../bld";   # List of where to look for the Perl modules
require "CAM_test.pm";


#
# Setup tests that will be run and create a new test object
#
my $cam = CAM_test->new;

#
# Process input arguments
#
$cam->process_args;

# BEGIN F_OPTIMIZATION_OVERRIDE HACK
# HACK:  F_OPTIMIZATION_OVERRIDE env var for dycore-dependent optimization
my $user_defined_fopt = "FALSE";
if ( defined($F_OPTIMIZATION_OVERRIDE) ) {
  $user_defined_fopt = "TRUE";
}
my $pgf90_opt_hack = "FALSE";
if (($cam->arch eq "Linux") && 
    (( ! defined($USER_FC)) || ($USER_FC eq "pgf90"))) {
  $pgf90_opt_hack = "TRUE";
  print "WARNING:  test-model.pl will reduce optimization levels for some\n";
  print "WARNING:  tests of the fv dycore due to problems with pgf90.\n";
}
# END F_OPTIMIZATION_OVERRIDE HACK

$cam->do_archivelog( "no" );  # Don't archive log files
my $platform = $cam->Platform;# Three letter description of platform
my @cprout;                   # Array of cprnc log files
my $DYNAMICS = "none"; 
$cam->CFG->setcfg( "DYNAMICS", $DYNAMICS );
#
# Set directories
#
my $ROOT_DIR  = "$SCRIPT_DIR";
$ROOT_DIR =~ s/\/models\/atm\/cam\/test\/system//;
my $MODEL_CPRDIR = "$ROOT_DIR/models/atm/cam/tools/cprnc"; # Comparision program
$cam->setenv( "MODEL_CPRDIR", $MODEL_CPRDIR );
$cam->CFG->setcfg( "CAMROOT", $ROOT_DIR );
#
# If environment variables set for specific configuration use it
#
if ( $PHYSICS =~ /./ )  { $cam->CFG->setcfg( "PHYSICS",  "$PHYSICS" ); }
if ( $OCEANMODEL =~ /./ ) { $cam->CFG->setcfg( "OCEANMODEL", "$OCEANMODEL" ); }
if ( $ICEMODEL =~ /./ ) { $cam->CFG->setcfg( "ICEMODEL",    "$ICEMODEL" ); }
# Keep track of initial setting of SPMD mode
my $SPMD_prev;
if ( defined($SPMD) ) { 
  $SPMD_prev = $SPMD;
} else {
  $SPMD_prev = undef;
}
# Keep track of initial setting of SMP mode
my $SMP_prev;
if ( defined($SMP) ) { 
  $SMP_prev = $SMP;
} else {
  $SMP_prev = undef;
}
#
# Definition of filenames that will be used for tests
#
my $filewildcard = "history1_instant*-??-??.nc";
my $filespec     = "history1_instant%y-%m-%d.nc";
my $avg_file     = "history2_avg.nc";
my $firstfile;
#
# Build support programs  (cprnc to compare history files)
#
$cam->setup_directories( "support" ); # setup the directories to build the support programs
my $CASE_DIR = $cam->env( "CASE_DIR" );
$cam->setenv( "EXEDIR", $CASE_DIR );    # Location to build cprnc and makdep
$cam->setenv( "VPATH", $MODEL_CPRDIR ); # Location of source for cprnc
$cam->clean( "$MODEL_CPRDIR" );       # Clean the source directory
my $BUILD_DIR = $cam->env( "BUILD_DIR" );
$cam->chdir( "$BUILD_DIR/cprncobj" ); # Change to build directory for cprnc
$cam->build_support( ("VPATH", "EXEDIR", "LIB_NETCDF", "INC_NETCDF") );
$cam->setenv( "EXEDIR", undef );
$cam->setenv( "VPATH", undef );
$cam->chdir( "$SCRIPT_DIR" ); # Change to build directory for cprnc
$cam->continuation_mark( "0_nothing_done" );  # Mark how far script has proceeded
#--------------------------------------------------------------------------
# Loop over each dynamics
#--------------------------------------------------------------------------
foreach my $dyn ( $cam->dynamics ) {
  my $trace_gas = ".true.";
  $cam->CFG->setcfg( "PERGRO", "FALSE" );  # Error growth mode is off
  $DYNAMICS     = $dyn; 
  $cam->CFG->setcfg( "DYNAMICS", $DYNAMICS );
  my $PLEV = $cam->config("PLEV");
  $cam->CFG->setcfg( "PLEV", $PLEV );
  my $RESOLUTION = $cam->config("RESOLUTION");
  $cam->CFG->setcfg( "RESOLUTION", $RESOLUTION );
  if ( $dyn eq "fv" ) {
    $cam->CFG->setcfg( "TWOD_YZ", $cam->config("TWOD_YZ") );
  } else {
    $cam->CFG->unsetcfg( "TWOD_YZ" );
  }
  my $adiabat = "false";
  &test_model_set_constits( $cam, $trace_gas );
  #
  # Variables that will be used later
  #
  my $MODEL_BLDDIR = undef;
  my $MODEL_EXEDIR = undef;
  my $RUNTYPE = undef;
  my $CASE = undef;
  $cam->continuation_mark( "0_nothing_done" );       # Mark how far script has proceeded
  #---------------------------------------------------------------------------
  # Configure/build the model with debugging on to catch errors
  #---------------------------------------------------------------------------
  my %debug_files;
  my $desc;
  my %save_exedir;
  foreach $desc ( $cam->tests("debug") ) {
    #
    # Setup configuration
    #
    $SPMD  = $cam->config("SPMD",$desc); # SPMD on or off
    $cam->CFG->setcfg( "SPMD", $SPMD );
    $SPMD  =~ /(.)/;  # Extract first character from SPMD
    my $SPMDC = $1;
    $SMP  = $cam->config("SMP",$desc); # SMP on or off
    $cam->CFG->setcfg( "SMP", $SMP );
    $SMP  =~ /(.)/;  # Extract first character from SMP
    my $SMPC = $1;
    $cam->CFG->setcfg( "DEBUG", "TRUE" ); # Debug (compiler error checking) turned on
    my $CASE  = "debugSPMD${SPMDC}_SMP${SMPC}${platform}${DYNAMICS}${RESOLUTION}L${PLEV}";  # Case name
    $cam->setenv( "CASE", $CASE );
    $cam->CFG->unsetcfg( "MODEL_BLDDIR" );  # Model build directory (set in setup_directories)
    $cam->CFG->unsetcfg( "MODEL_EXEDIR" );  # Model execution directory (set in setup_directories)
    #
    # Change directory, setup build/run directories, build model
    #
    $cam->chdir( $SCRIPT_DIR ); 
    $cam->setup_directories; # Setup the build/run/source directories
    $save_exedir{$desc} = $cam->CFG->cfg( "MODEL_EXEDIR" );
    $cam->chdir( "MODEL_BLDDIR" );# Change to build directory
    $cam->configure;               # Run the configure script
    $cam->make;                    # Build the model executable
    #---------------------------------------------------------------------------
    # Run the model with debugging on
    #---------------------------------------------------------------------------
    $cam->chdir( "MODEL_EXEDIR" ); # Change to exedir
    #
    # Define the namelist
    #
    &test_model_reset_namelists( $cam, $desc, $adiabat, $filespec, $avg_file );
    if ( ($RESOLUTION eq "32x64") || ($RESOLUTION eq "4x5") || ($RESOLUTION eq "8x16") ||
         ($RESOLUTION eq "64x128") ) {
      $cam->setenv( "NCDATA_VERS", 2 );
    }
    $main::CAMEXP{'nestep'}    = 3;
    #
    # Build namelist, run model, save files, and clean afterwards
    #
    $RUNTYPE = "initial";
    $cam->setenv( "RUNTYPE", $RUNTYPE );
    $cam->run_time_env;            # Run-time env vars
    $cam->build_namelist( 1 );     # Build the model namelist
    my $log_file = $cam->env( "LOG_DIR" ) . "/$desc.log";
    $cam->run( "$desc: Run $RUNTYPE simulation with DEBUG on for # steps=: ".
               $main::CAMEXP{'nestep'}, $log_file );
    $MODEL_EXEDIR = $cam->CFG->cfg( "MODEL_EXEDIR" );
    ($debug_files{$desc}) = glob("$MODEL_EXEDIR/$filewildcard");
    if ( ! $cam->do_continue && ! defined($debug_files{$desc}) ) {
      $cam->die( "Output file does not exist for wildcard: ". $filewildcard );
    }
    if ( $cam->do_clean ) { $cam->clean; }   # Clean out BLD and EXE directories
    #
    # Compare the history files with and without SPMD-mode
    #
    if ( $desc eq '03_initial_debug_run_nonSPMD'
      && -f $debug_files{'01_initial_debug_run_SPMD'} 
      && -f $debug_files{'03_initial_debug_run_nonSPMD'} ) {
      $cam->compare_files( $debug_files{'01_initial_debug_run_SPMD'}, 
                         $debug_files{'03_initial_debug_run_nonSPMD'}, 
                         "$MODEL_EXEDIR/${desc}_1to3.cprout", 
                         "SPMD and non-SPMD mode are are not bit-for-bit" );
    }
    #
    # Compare the history files with and without SMP-mode
    #
    if ( $platform !~ /aix/ ) {
      if ( $desc eq '02_initial_debug_run_SPMDnoSMD'
        && -f $debug_files{'01_initial_debug_run_SPMD'} 
        && -f $debug_files{'02_initial_debug_run_SPMDnoSMD'} ) {
        $cam->compare_files( $debug_files{'01_initial_debug_run_SPMD'}, 
                             $debug_files{'02_initial_debug_run_SPMDnoSMD'}, 
                             "$MODEL_EXEDIR/${desc}_1to2.cprout", 
                             "SMP and non-SMP mode are are not bit-for-bit" );
      }
    }
    $cam->continuation_mark( $desc );  # Mark how far script has proceeded
  }
  #
  # Clean all debug datafiles out now
  #
  if ( $cam->do_clean ) { 
    foreach $desc ( $cam->tests("debug") ) {
      if ( defined($save_exedir{$desc}) ) {
        $cam->CFG->setcfg( "MODEL_EXEDIR", $save_exedir{$desc} );
        $cam->clean( "alldata" );
       }
    }
  }
  $SPMD = $SPMD_prev;
  if ( defined($SPMD) ) {
    $cam->CFG->setcfg( "SPMD", $SPMD );
  } else {
    $cam->CFG->unsetcfg( "SPMD" );
  }
  $SMP = $SMP_prev;
  if ( defined($SMP) ) {
    $cam->CFG->setcfg( "SMP", $SMP );
  } else {
    $cam->CFG->unsetcfg( "SMP" );
  }
  #---------------------------------------------------------------------------
  # Now configure/build/run with full physics to test restarts different PE's etc.
  #---------------------------------------------------------------------------
  my %test_files;
  if ( $cam->tests("restart") ) {
    #
    # Setup configuration
    #
    $cam->CFG->setcfg( "DEBUG", "FALSE" ); # Debug (compiler error checking) turned off
    # Use a long casename to ensure that long casenames function correctly
    $CASE = "test$platform${DYNAMICS}${RESOLUTION}L${PLEV}_0123456789012";  # Case name
    $cam->setenv( "CASE", $CASE );
    $cam->CFG->unsetcfg( "MODEL_BLDDIR" );  # Model build directory (set in setup_directories)
    $cam->CFG->unsetcfg( "MODEL_EXEDIR" );  # Model execution directory (set in setup_directories)
    #
    # Change directory, setup build/run directories, build model
    #
    $cam->chdir( $SCRIPT_DIR ); 
    $cam->setup_directories; # Setup the build/run/source directories
    $save_exedir{'restart'} = $cam->CFG->cfg( "MODEL_EXEDIR" );
    $cam->chdir( "MODEL_BLDDIR" ); # Change to build directory
    $cam->configure;               # Run the configure script
    # BEGIN F_OPTIMIZATION_OVERRIDE HACK
    if ( $dyn eq "fv" ) {
      if ( $user_defined_fopt eq "FALSE" ) {
        if ( $pgf90_opt_hack eq "TRUE" ) {
          $F_OPTIMIZATION_OVERRIDE = -O0;
          print "WARNING:  test-model.pl reduced pgf90 optimization to -O0\n";
        }
      }
    }
    # END F_OPTIMIZATION_OVERRIDE HACK
    $cam->make;                    # Build the model executable
    # BEGIN F_OPTIMIZATION_OVERRIDE HACK
    if ( $dyn eq "fv" ) {
      if ( $user_defined_fopt eq "FALSE" ) {
        if ( $pgf90_opt_hack eq "TRUE" ) {
          undef $F_OPTIMIZATION_OVERRIDE;
        }
      }
    }
    # END F_OPTIMIZATION_OVERRIDE HACK
    #---------------------------------------------------------------------------
    # Run the model configuration
    # (ensure model starts, restarts and restarts are bit-for-bit)
    # (change PE's on restart to make sure restarts not PE dependent)
    #---------------------------------------------------------------------------
    $cam->chdir( "MODEL_EXEDIR" ); # Change to exedir
  }
  foreach $desc ( $cam->tests("restart") ) {
    #
    # Define the namelist
    #
    $RUNTYPE = $cam->config("RUNTYPE", $desc );
    $cam->setenv( "RUNTYPE", $RUNTYPE );
    &test_model_reset_namelists( $cam, $desc, $adiabat, $filespec, $avg_file );
    $main::CAMEXP{'nhtfrq(1)'} = $cam->config("nestep", '04_initial' );
    $main::CAMEXP{'nestep'}    = $cam->config("nestep", $desc );
    #
    # Build namelist, run model, save files
    #
    if ( $desc eq "05_restart" ) { $cam->lower_PEs; }
    $cam->run_time_env;            # Run-time env vars
    $cam->build_namelist( 1 );     # Build the model namelist
    if ( -f glob("$filewildcard") && ($RUNTYPE ne "restart") ) { 
      $cam->exec( "/bin/rm $filewildcard" ); # Remove file before run
    }
    if ( -f $test_files{$desc}{h} ) {
      $cam->exec( "/bin/rm $test_files{$desc}{h}" );  # Remove file before run
    }
    if ( -f $test_files{$desc}{avg} ) {
      $cam->exec( "/bin/rm $test_files{$desc}{avg}" ); # Remove file before run
    }
    my $log_file = $cam->env( "LOG_DIR" ) . "/$desc.log"; 
    $cam->run( "$desc: Run $RUNTYPE simulation for # steps=: ".$cam->config("nestep", $desc ), 
               "$log_file" );
    # Compare primary files
    $MODEL_EXEDIR = $cam->CFG->cfg( "MODEL_EXEDIR" );
    $test_files{$desc}{h} = "$MODEL_EXEDIR/history1.$desc.nc";
    ($firstfile) = glob("$filewildcard");
    if ( ! $cam->do_continue && ! defined($firstfile) ) {
      $cam->die( "Output file does not exist for wildcard: ". $filewildcard );
    }
    $cam->exec( "/bin/cp $firstfile $test_files{$desc}{h}" );
    if ( $desc eq "06_initial_compare_to_restart" ) { 
      $cam->compare_files( $test_files{$desc}{h}, $test_files{'05_restart'}{h}, 
         "$MODEL_EXEDIR/${desc}_first.cprout", "Restart or change in PE not exact" );
    }
    # Compare Secondary avg file
    $test_files{$desc}{avg} = "$MODEL_EXEDIR/history2_avg.$desc.nc";
    $cam->exec( "/bin/cp $avg_file $test_files{$desc}{avg}" );
    if ( $desc eq "06_initial_compare_to_restart" ) { 
      $cam->compare_files( $test_files{$desc}{avg}, $test_files{'05_restart'}{avg}, 
         "$MODEL_EXEDIR/${desc}_avg.cprout", "Restart or change in PE not exact for secondary files" );
      $cam->save_timing( $desc, $log_file );
    }
    if ( $desc eq "05_restart" ) { $cam->reset_PEs; }
    $cam->continuation_mark( $desc );  # Mark how far script has proceeded
  }
  if ( $cam->do_clean ) { $cam->clean; }   # Clean out BLD and EXE directories
  #---------------------------------------------------------------------------
  # Now configure/build/run the code for special physics configurations (i.e. SOM)
  #---------------------------------------------------------------------------
  my $adiabat = "false";
  my %files;
  foreach $desc ( $cam->tests("physics") ) {
    #
    # Setup configuration
    #
    $cam->CFG->setcfg( "OCEANMODEL", $cam->config("OCEANMODEL",$desc) );
    $cam->CFG->setcfg( "DEBUG", "TRUE" ); # Debug (compiler error checking) turned on
    $CASE = "${desc}${platform}${DYNAMICS}${RESOLUTION}L${PLEV}";  # Case name
    $cam->setenv( "CASE", $CASE );
    $cam->CFG->unsetcfg( "MODEL_BLDDIR" );  # Model build directory (set in setup_directories)
    $cam->CFG->unsetcfg( "MODEL_EXEDIR" );  # Model execution directory (set in setup_directories)
    #
    # Change directory, setup build/run directories, build model
    #
    $cam->chdir( $SCRIPT_DIR ); 
    $cam->setup_directories( ); # Setup the build/run/source directories
    $cam->chdir( "MODEL_BLDDIR" );   # Change to build directory
    $MODEL_EXEDIR = $cam->CFG->cfg( "MODEL_EXEDIR" );
    $save_exedir{$desc} = $cam->CFG->cfg( "MODEL_EXEDIR" );
    $files{$adiabat}{$desc} = "$MODEL_EXEDIR/history1.$desc.nc";
    $cam->configure;               # Run the configure script
    $cam->make;                    # Build the model executable
    #---------------------------------------------------------------------------
    # Run the model
    #---------------------------------------------------------------------------
    $cam->chdir( "MODEL_EXEDIR" ); # Change to exedir
    #
    # Define the namelist
    #
    &test_model_reset_namelists( $cam, $desc, $adiabat, $filespec, $avg_file );
    $main::CAMEXP{'nestep'}    = $cam->config("nestep", '04_initial' );
    #
    # Build namelist, run model, save files, and clean afterwards
    #
    $RUNTYPE = "initial";
    $cam->setenv( "RUNTYPE", $RUNTYPE );
    $cam->run_time_env;            # Run-time env vars
    $cam->build_namelist( 1 );     # Build the model namelist
    my $log_file = $cam->env( "LOG_DIR" ) . "/$desc.log";
    $cam->run( "$desc: Run $RUNTYPE simulation with $desc for # steps=: ".
               $main::CAMEXP{'nestep'}, $log_file );
    ($firstfile) = glob("$filewildcard");
    if ( ! $cam->do_continue && ! defined($firstfile) ) {
      $cam->die( "Output file does not exist for wildcard: ". $filewildcard );
    }
    $cam->exec( "/bin/cp $firstfile $files{$adiabat}{$desc}" );
    if ( $cam->do_clean ) { $cam->clean; }   # Clean out BLD and EXE directories
    $cam->continuation_mark( $desc );  # Mark how far script has proceeded
    #
    # Reset configuration
    #
    if ( $OCEANMODEL =~ /./ ) { $cam->CFG->setcfg( "OCEANMODEL", "$OCEANMODEL" ); }
    else { $cam->CFG->unsetcfg( "OCEANMODEL" ); }
  }
  #---------------------------------------------------------------------------
  # Now configure/build/run the control code without error-growth and with full-physics on
  #---------------------------------------------------------------------------
  foreach $desc ( $cam->tests("control") ) {
    #
    # Setup configuration
    #
    $cam->CFG->setcfg( "DEBUG", "TRUE" ); # Debug (compiler error checking) turned on
    $cam->CFG->setcfg( "OCEANMODEL", $cam->config("OCEANMODEL",$desc) );
    $CASE = "control${platform}${DYNAMICS}${RESOLUTION}L${PLEV}";  # Case name
    $cam->setenv( "CASE", $CASE );
    $cam->CFG->unsetcfg( "MODEL_BLDDIR" );  # Model build directory (set in setup_directories)
    $cam->CFG->unsetcfg( "MODEL_EXEDIR" );  # Model execution directory (set in setup_directories)
    #
    # Change directory, setup build/run directories, build model
    #
    $cam->chdir( $SCRIPT_DIR ); 
    $cam->setup_directories( "compare" ); # Setup the build/run/source directories
    $cam->chdir( "MODEL_BLDDIR" );   # Change to build directory
    $MODEL_EXEDIR = $cam->CFG->cfg( "MODEL_EXEDIR" );
    $files{$adiabat}{$desc} = "$MODEL_EXEDIR/history1.$desc.nc";
    $cam->configure;               # Run the configure script
    $cam->make;                    # Build the model executable
    #---------------------------------------------------------------------------
    # Run the control case model
    #---------------------------------------------------------------------------
    $cam->chdir( "MODEL_EXEDIR" ); # Change to exedir
    #
    # Define the namelist
    #
    &test_model_reset_namelists( $cam, $desc, $adiabat, $filespec, $avg_file );
    $main::CAMEXP{'nestep'}    = $cam->config("nestep", '06_initial_compare_to_restart' );
    #
    # Build namelist, run model, save files, and clean afterwards
    #
    $RUNTYPE = "initial";
    $cam->setenv( "RUNTYPE", $RUNTYPE );
    $cam->run_time_env;            # Run-time env vars
    $cam->build_namelist( 1 );     # Build the model namelist
    my $log_file = $cam->env( "LOG_DIR" ) . "/$desc.log";
    $cam->run( "$desc: Run $RUNTYPE simulation with control library for # steps=: ".
               $main::CAMEXP{'nestep'}, $log_file );
    ($firstfile) = glob("$filewildcard");
    if ( ! $cam->do_continue && ! defined($firstfile) ) {
      $cam->die( "Output file does not exist for wildcard: ". $filewildcard );
    }
    $cam->exec( "/bin/cp $firstfile $files{$adiabat}{$desc}" );
    if ( $cam->do_clean ) { $cam->clean; }   # Clean out BLD and EXE directories
    $cam->save_timing( $desc, $log_file );
    $cam->continuation_mark( $desc );  # Mark how far script has proceeded
    #
    # Reset configuration
    #
    if ( $OCEANMODEL =~ /./ ) { $cam->CFG->setcfg( "OCEANMODEL", "$OCEANMODEL" ); }
    else { $cam->CFG->unsetcfg( "OCEANMODEL" ); }
  }
  #
  # Compare files from control cases to experiment cases
  #
  my $MODEL_EXEDIR = $cam->CFG->cfg( "MODEL_EXEDIR" );
  if ( -f $test_files{'06_initial_compare_to_restart'}{'h'} 
    && -f $files{$adiabat}{'08_control_nonpert'} ) {
    push @cprout, "$MODEL_EXEDIR/con.$desc.cprout";
    $cam->compare_files( $test_files{'06_initial_compare_to_restart'}{'h'}, 
         $files{$adiabat}{'08_control_nonpert'}, $cprout[$#cprout] );
  }
  if ( -f $test_files{'07_SOM'}{'h'} && -f $files{$adiabat}{'09_control_SOM'} ) {
    push @cprout, "$MODEL_EXEDIR/con.$desc.cprout";
    $cam->compare_files( $test_files{'07_SOM'}{'h'}, 
         $files{$adiabat}{'09_control_SOM'}, $cprout[$#cprout] );
  }
  #
  # Clean data files out
  #
  if ( $cam->do_clean ) {
    if ( defined($save_exedir{'restart'}) ) {
      $cam->CFG->setcfg( "MODEL_EXEDIR", $save_exedir{'restart'} );
      $cam->clean( "alldata" );
    }
    foreach $desc ( $cam->tests("physics") ) {
      if ( defined($save_exedir{$desc}) ) {
        $cam->CFG->setcfg( "MODEL_EXEDIR", $save_exedir{$desc} );
        $cam->clean( "alldata" );
      }
    }
    foreach $desc ( $cam->tests("control") ) {
      if ( defined($save_exedir{$desc}) ) {
        $cam->CFG->setcfg( "MODEL_EXEDIR", $save_exedir{$desc} );
        $cam->clean( "alldata" );
      }
    }
  }
  #---------------------------------------------------------------------------
  # Now configure/build/run an error growth test
  # (Run both with adiabatic on and off)
  #---------------------------------------------------------------------------
  #
  # Turn trace_gas off, and read-trace on, to check a slightly different configuration
  #
  my $roundoff = $cam->config( "roundoff" );
  $trace_gas = ".false.";
  if ( $cam->tests("errgro") ) {
    #
    # Setup configuration
    #
    &test_model_set_constits( $cam, $trace_gas );
    $cam->CFG->setcfg( "PERGRO", "TRUE" );  # Error growth on
    $cam->CFG->setcfg( "DEBUG", "FALSE" );  # Debug (compiler-checking) off
    $CASE = "errgro$platform${DYNAMICS}${RESOLUTION}L${PLEV}";# Case name
    $cam->setenv( "CASE", $CASE );
    $cam->CFG->unsetcfg( "MODEL_BLDDIR" );  # Model build directory (set in setup_directories)
    $cam->CFG->unsetcfg( "MODEL_EXEDIR" );  # Model execution directory (set in setup_directories)
    #
    # Change directory, setup build/run directories, build model
    #
    $cam->chdir( $SCRIPT_DIR ); 
    $cam->setup_directories; # Setup the build/run/source directories
    $save_exedir{'errgro'} = $cam->CFG->cfg( "MODEL_EXEDIR" );
    $cam->chdir( "MODEL_BLDDIR" );   # Change to build directory
    $cam->configure;                 # Run the configure script
    $cam->run_time_env;              # Run-time env vars
    if ( $dyn eq "fv" ) {
      if ( $user_defined_fopt eq "FALSE" ) {
        if ( $pgf90_opt_hack eq "TRUE" ) {
          $F_OPTIMIZATION_OVERRIDE = -O1;
          print "WARNING:  test-model.pl reduced pgf90 optimization to -O1\n";
        }
      }
    }
    $cam->make;                    # Build the model executable
    if ( $dyn eq "fv" ) {
      if ( $user_defined_fopt eq "FALSE" ) {
        if ( $pgf90_opt_hack eq "TRUE" ) {
          undef $F_OPTIMIZATION_OVERRIDE;
        }
      }
    }
    #
    # Run the model
    #
    $cam->chdir( "MODEL_EXEDIR" );   # Change to exedir
  }
  foreach $desc ( $cam->tests("errgro") ) {
    my $adiabat = $cam->config("adiabatic",$desc );
    my $pertlim = $cam->config("pertlim", $desc );
    if ( ! defined($adiabat) || ! defined($pertlim) ) { 
      die "Error growth: $0 script error adiabat or pertlim not set"; 
    }
    #
    # Define namelist
    #
    &test_model_reset_namelists( $cam, $desc, $adiabat, $filespec, $avg_file );
    $main::CAMEXP{'pertlim'} = $pertlim;
    #
    # Build namelist, run model, save files
    #
    $RUNTYPE = "initial";
    $cam->setenv( "RUNTYPE", $RUNTYPE );
    $cam->run_time_env;            # Run-time env vars
    $cam->build_namelist( 1 );     # Build the model namelist
    $MODEL_EXEDIR = $cam->CFG->cfg( "MODEL_EXEDIR" );
    $files{$adiabat}{$pertlim} = "$MODEL_EXEDIR/history1.$desc.nc";
    if ( -f glob("$filewildcard") ) { 
      $cam->exec( "/bin/rm $filewildcard" ); # Remove file before run
    }
    if ( -f $files{$adiabat}{$pertlim} ) {
      $cam->exec( "/bin/rm $files{$adiabat}{$pertlim}" ); # Remove file before run
    }
    my $log_file = $cam->env( "LOG_DIR" ) . "/$desc.log";
    $cam->run( "$desc: Run error growth adiabatic set to $adiabat with perturbation set to: " 
               . "$pertlim", $log_file );
    ($firstfile) = glob("$filewildcard");
    if ( ! $cam->do_continue && ! defined($firstfile) ) {
      $cam->die( "Output file does not exist for wildcard: ". $filewildcard );
    } else {
      $cam->exec( "/bin/cp $firstfile $files{$adiabat}{$pertlim}" );
    }
    #
    # Compare files
    #
    if ( $pertlim == $roundoff ) {
      push @cprout, "$MODEL_EXEDIR/$CASE.$desc.cprout";
      if ( $cam->do_remote_errg ) {
        $cam->compare_files( $cam->remote_baseline_errg( $files{$adiabat}{0} ), 
                             $files{$adiabat}{$roundoff}, $cprout[$#cprout] );
      } else {
        $cam->compare_files( $files{$adiabat}{0}, $files{$adiabat}{$roundoff}, 
                             $cprout[$#cprout] );
      }
    } else {
      $cam->save_timing( $desc, $log_file );
    }
    $cam->continuation_mark( "$desc" );  # Mark how far script has proceeded
  }
  #---------------------------------------------------------------------------
  # Compare to control source directory.
  # always compare if files exist, only rerun if they don't.
  #---------------------------------------------------------------------------
  #
  # Setup configuration
  #
  $trace_gas = ".false.";
  &test_model_set_constits( $cam, $trace_gas );
  $cam->CFG->setcfg( "PERGRO", "TRUE" );  # Error growth on
  $cam->CFG->setcfg( "DEBUG", "FALSE" ); # Debug (compiler error checking) turned off
  $CASE   = "controlerrg$platform${DYNAMICS}${RESOLUTION}L${PLEV}";  # Case name
  $cam->setenv( "CASE", $CASE );
  $cam->CFG->unsetcfg( "MODEL_BLDDIR" );  # Model build directory (set in setup_directories)
  $cam->CFG->unsetcfg( "MODEL_EXEDIR" );  # Model execution directory (set in setup_directories)
  #
  # Change directory, setup build/run directories, build model
  #
  $cam->chdir( $SCRIPT_DIR ); 
  $cam->setup_directories( "compare" );# Setup the directories
  if ( $cam->tests("conterrgro") ) {
    $cam->chdir( "MODEL_BLDDIR" );   # Change to build directory
    $cam->configure;                 # Run the configure script
    $cam->make;                      # Build the model executable
  }
  #
  # Run the model
  #
  $cam->chdir( "MODEL_EXEDIR" );     # Change to exedir
  foreach $desc ( $cam->tests("conterrgro") ) {
    my $adiabat = $cam->config("adiabatic",$desc );
    my $pertlim = $cam->config("pertlim", $desc );
    if ( ! defined($adiabat) || ! defined($pertlim) ) { 
      die "Control: $0 script error adiabat or pertlim not set"; 
    }
    #
    # Define the namelist
    #
    &test_model_reset_namelists( $cam, $desc, $adiabat, $filespec, $avg_file );
    $main::CAMEXP{'pertlim'} = $pertlim;
    #
    # Build namelist, run model, save files, and clean afterwards
    #
    $MODEL_EXEDIR = $cam->CFG->cfg( "MODEL_EXEDIR" );
    $files{$adiabat}{"control$pertlim"} = "$MODEL_EXEDIR/history1.$desc.nc";
    $cam->run_time_env;            # Run-time env vars
    $cam->build_namelist( 1 );          # Build the model namelist
    if ( -f glob("$filewildcard") ) { 
      $cam->exec( "/bin/rm $filewildcard" ); # Remove file before run
    }
    if ( -f $files{$adiabat}{"control$pertlim"} ) { 
      $cam->exec( "/bin/rm ".$files{$adiabat}{"control$pertlim"} ); # Remove file before run
    }
    my $log_file = $cam->env( "LOG_DIR" ) . "/$desc.log";
    $cam->run( "$desc: Run control case to compare to, adiabatic set to :".$adiabat.
               " Pertlim: $pertlim ", "$log_file" );
    ($firstfile) = glob("$filewildcard");
    if ( ! $cam->do_continue && ! defined($firstfile) ) {
      $cam->die( "Output file does not exist for wildcard: ". $filewildcard );
    }
    $cam->exec( "/bin/cp $firstfile ".$files{$adiabat}{"control$pertlim"} );
    #
    # Compare files
    #
    push( @cprout, "$MODEL_EXEDIR/$CASE.$desc.cprout" );
    if ( $pertlim == $roundoff ) {
      if ( -f $files{$adiabat}{"control$roundoff"} && -f $files{$adiabat}{control0} ) {
        $cam->compare_files( $files{$adiabat}{control0}, $files{$adiabat}{"control$roundoff"}, $cprout[$#cprout] );
      }
    } else {
      if ( -f $files{$adiabat}{control0} && -f $files{$adiabat}{0} ) {
        $cam->compare_files( $files{$adiabat}{0}, $files{$adiabat}{control0}, $cprout[$#cprout] );
      }
    }
    $cam->save_timing( $desc, $log_file );
    $cam->continuation_mark( "$desc" );  # Mark how far script has proceeded
  }
  #
  # Clean out data files
  #
  if ( $cam->do_clean ) {
    if ( $cam->tests("conterrgro") ) {
      $cam->clean;
      $cam->clean( "alldata" );
    }
    if ( defined($save_exedir{'errgro'}) ) {
      $cam->CFG->setcfg( "MODEL_EXEDIR", $save_exedir{'errgro'} );
      $cam->clean( "alldata" );
    }
  }
  #---------------------------------------------------------------------------
  # Configure/build the model for special configurations
  #---------------------------------------------------------------------------
  my %debug_files;
  my $desc;
  $cam->CFG->setcfg( "PERGRO", "FALSE" );  # Error growth mode off
  if ( $USER_FC eq "lf95" ) {
    $cam->CFG->setcfg( "DEBUG", "FALSE" ); # Debug (compiler error checking) turned off
  } else {
    $cam->CFG->setcfg( "DEBUG", "TRUE" ); # Debug (compiler error checking) turned on
  }
  my $PHYSICSold   = $PHYSICS;
  my $OCEANMODELold = $OCEANMODEL;
  my $LANDMODELold = $LANDMODEL;
  my $ICEMODELold  = $ICEMODEL;
  foreach $desc ( $cam->tests("special") ) {
    if ( $DYNAMICS eq $cam->config("DYNAMICS", $desc ) ) {
      $LANDMODEL = $cam->config("LANDMODEL", $desc );
      $PHYSICS = $cam->config("PHYSICS", $desc );
      $OCEANMODEL = $cam->config("OCEANMODEL", $desc );
      $ICEMODEL = $cam->config("ICEMODEL", $desc );
      if ( $LANDMODEL ne $LANDMODELold || $PHYSICS ne $PHYSICSold || $OCEANMODEL ne $OCEANMODELold
      ||   $ICEMODEL ne $ICEMODELold ) {
        #
        # Setup configuration
        #
        if ( defined($LANDMODEL) ) {
          $cam->CFG->setcfg( "LANDMODEL", $LANDMODEL );
        } else {
          $cam->CFG->unsetcfg( "LANDMODEL" );
        }
        if ( defined($ICEMODEL) ) {
          $cam->CFG->setcfg( "ICEMODEL", $ICEMODEL );
        } else {
          $cam->CFG->unsetcfg( "ICEMODEL" );
        }
        if ( defined($PHYSICS) ) {
          $cam->CFG->setcfg( "PHYSICS", $PHYSICS );
        } else {
          $cam->CFG->unsetcfg( "PHYSICS" );
        }
        if ( defined($OCEANMODEL) ) {
          $cam->CFG->setcfg( "OCEANMODEL", $OCEANMODEL );
        } else {
          $cam->CFG->unsetcfg( "OCEANMODEL" );
        }
        $RESOLUTION = $cam->config("RESOLUTION", $desc );
        if ( ! defined($RESOLUTION) ) { $RESOLUTION = $cam->config("RESOLUTION"); }
        $cam->CFG->setcfg( "RESOLUTION", $RESOLUTION );
        $PLEV = $cam->config("PLEV", $desc );
        if ( ! defined($PLEV) ) {   $PLEV = $cam->config("PLEV"); }
        $cam->CFG->setcfg( "PLEV", $PLEV );
        if ( defined($PHYSICS) && ($PHYSICS eq "ccm366") ) {
          $trace_gas = ".false.";
        } else {
          $trace_gas = ".true.";
        }
        &test_model_set_constits( $cam, $trace_gas );
        $CASE  = "$desc${platform}${DYNAMICS}${RESOLUTION}L${PLEV}";  # Case name
        $cam->setenv( "CASE", $CASE );
        $cam->CFG->unsetcfg( "MODEL_BLDDIR" );  # Model build directory (set in setup_directories)
        $cam->CFG->unsetcfg( "MODEL_EXEDIR" );  # Model execution directory (set in setup_directories)
        #
        # Change directory, setup build/run directories, build model
        #
        $cam->chdir( $SCRIPT_DIR ); 
        $cam->setup_directories; # Setup the build/run/source directories
        $cam->chdir( "MODEL_BLDDIR" ); # Change to build directory
        $cam->configure;               # Run the configure script
        $cam->run_time_env;            # Run-time env vars
        $cam->make;                    # Build the model executable
        #---------------------------------------------------------------------------
        # Run the model with debugging on
        #---------------------------------------------------------------------------
        $cam->chdir( "MODEL_EXEDIR" ); # Change to exedir
        #
        # Define namelist
        #
        &test_model_reset_namelists( $cam, $desc, $adiabat, $filespec, $avg_file );
        $main::CAMEXP{'nestep'} = 3;
        #
        # Build namelist, run model, save files, and clean afterwards
        #
        $RUNTYPE = "initial";
        $cam->setenv( "RUNTYPE", $RUNTYPE );
        $cam->run_time_env;            # Run-time env vars
        $cam->build_namelist( 1 );     # Build the model namelist
        my $log_file = $cam->env( "LOG_DIR" ) . "/$desc.log";
        $cam->run( 
        "$desc: Run $RUNTYPE simulation: $PHYSICS physics,$LANDMODEL, $ICEMODEL and $OCEANMODEL: DEBUG on for # steps=: ".
                   $main::CAMEXP{'nestep'}, "$log_file" );
        $MODEL_EXEDIR = $cam->CFG->cfg( "MODEL_EXEDIR" );
        ($debug_files{$desc}) = glob("$MODEL_EXEDIR/$filewildcard");
        if ( ! $cam->do_continue && ! defined($debug_files{$desc}) ) {
          $cam->die( "Output file does not exist for wildcard: ". $filewildcard );
        }
        if ( $cam->do_clean ) { 
          $cam->clean;     # Clean out BLD and EXE directories
          $cam->clean( "alldata" ); 
        }
        $cam->continuation_mark( $desc );  # Mark how far script has proceeded
        #
        # Set back to original values, set internal values in case
        # you are setting back to the default undefined values.
        #
        $LANDMODEL = $LANDMODELold;
        $ICEMODEL  = $ICEMODELold;
        $PHYSICS   = $PHYSICSold;
        $OCEANMODEL = $OCEANMODELold;
        if ( defined($LANDMODEL) ) {
          $cam->CFG->setcfg( "LANDMODEL", $LANDMODEL );
        } else {
          $cam->CFG->unsetcfg( "LANDMODEL" );
        }
        if ( defined($ICEMODEL) ) {
          $cam->CFG->setcfg( "ICEMODEL", $ICEMODEL );
        } else {
          $cam->CFG->unsetcfg( "ICEMODEL" );
        }
        if ( defined($PHYSICS) ) {
          $cam->CFG->setcfg( "PHYSICS", $PHYSICS );
        } else {
          $cam->CFG->unsetcfg( "PHYSICS" );
        }
        if ( defined($OCEANMODEL) ) {
          $cam->CFG->setcfg( "OCEANMODEL", $OCEANMODEL );
        } else {
          $cam->CFG->unsetcfg( "OCEANMODEL" );
        }
      }
    }
  }
}
print "Error growth comparision files are in: @cprout\n";
$cam->report_compare( @cprout );    # Report on success of comparisions
$cam->continuation_mark( "done" );  # Mark that script is completed


sub test_model_reset_namelists {
#
# Reset the namelists (this makes sure namelist is what it should be each time)
#
  my $cam = shift;
  my $desc = shift;
  my $adiabat = shift;
  my $filespec = shift;
  my $avg_file = shift;

  #
  # Undefine the namelists
  #
  %main::CAMEXP  = {};
  %main::LSMEXP  = {};
  %main::CLMEXP  = {};
  %main::MPRUN2D = {};
  %main::MPRUN2D = {}; %main::LSMEXP  = {}; %main::CLMEXP = {}; # Repeat so "use diagnostics" won't complain.
  foreach my $key ( keys(%main::CAMEXP) ) {
    $main::CAMEXP{$key} = undef;
  }
  foreach my $key ( keys(%main::LSMEXP) ) {
    $main::LSMEXP{$key} = undef;
  }
  foreach my $key ( keys(%main::CLMEXP) ) {
    $main::CLMEXP{$key} = undef;
  }
  foreach my $key ( keys(%main::MPRUN2D) ) {
    $main::MPRUN2D{$key} = undef;
  }
  #
  # Set values used by most
  #
  my $DYNAMICS = $cam->CFG->cfg("DYNAMICS");
  $main::CAMEXP{'readtrace'}    = ".false.";
  my $PHYSICS = $cam->CFG->cfg("PHYSICS");
  if ( $PHYSICS ne "ccm366" ) { 
    $main::CAMEXP{'trace_gas'}  = ".true.";
  } else {
    $main::CAMEXP{'trace_gas'}  = ".false.";
  }
  $main::CAMEXP{'adiabatic'}    = ".$adiabat.";
  $main::CAMEXP{'ctitle'}       = "\'$desc\'";
  $main::CAMEXP{'mss_irt'}      = 0;
  $main::CAMEXP{'ndens(1)'}     = 1;
  $main::CAMEXP{'ndens(2)'}     = 1;
  $main::CAMEXP{'nhtfrq(1)'}    = 1;
  $main::CAMEXP{'nhtfrq(2)'}    = $cam->config("nestep", "05_restart" );
  $main::CAMEXP{'hfilename_spec(1)'} = "\'$filespec\'";
  $main::CAMEXP{'hfilename_spec(2)'} = "\'$avg_file\'";
  #
  # Set NCDATA_VERS to the default
  #
  $cam->setenv( "NCDATA_VERS", undef );
  #
  # Keep restart pointer files out of $HOME
  #
  my $CASE_DIR = $cam->env( "CASE_DIR" );
  my $CASE = $cam->env( "CASE" );
  $main::CAMEXP{'rest_pfile'} = "\'$CASE_DIR/cam2.$CASE.rpointer\'";
  $main::CLMEXP{'rpntpath'}   = "\'$CASE_DIR/lnd.$CASE.rpointer\'";
  #
  # Make sure T and PS on tapes as needed for error growth
  #
  $main::CAMEXP{'fincl1(1)'}    = "\'T:I\'";
  $main::CAMEXP{'fincl1(2)'}    = "\'PS:I\'";
  #
  # Make sure secondary tapes are giving average values (rather than instant) so can 
  # verify continuing a average history tape mid-stream works correctly.
  #
  $main::CAMEXP{'fincl2(1)'}    = "\'T:A\'";
  $main::CAMEXP{'fincl2(2)'}    = "\'PS:A\'";
  $main::CAMEXP{'fincl2(3)'}    = "\'PSL:A\'";
  $main::CAMEXP{'fincl2(4)'}    = "\'U200:A\'";
  $main::CAMEXP{'fincl2(5)'}    = "\'U850:A\'";
  $main::CAMEXP{'fincl2(6)'}    = "\'V200:A\'";
  $main::CAMEXP{'fincl2(7)'}    = "\'V850:A\'";
  $main::CAMEXP{'fincl2(8)'}    = "\'T300:A\'";
  $main::CAMEXP{'fincl2(9)'}    = "\'T850:A\'";
  $main::CAMEXP{'fincl2(10)'}   = "\'Z300:A\'";
  $main::CAMEXP{'fincl2(11)'}   = "\'Z500:A\'";
  $main::CAMEXP{'fincl2(12)'}   = "\'Z700:A\'";
  #
  # Change date to just before a mid-month dataset update is necissary
  #
  $main::CAMEXP{'start_ymd'}    = 914;
  $main::CAMEXP{'start_tod'}    = 3600*24 - 1200;
  #
  # FV 2D parallel decomposition
  #
  if ( ($DYNAMICS eq "fv") && ($cam->CFG->cfg("TWOD_YZ") =~ /TRUE/) ) {
    if ( $SPMD_NODES !~ /./ ) {
      $cam->die( "SPMD_NODES not set when using FV 2D decomposition. ".
                 "Set the environment variable SPMD_NODES to the value of SPMD ".
                 "tasks to use for this configuration" );
    }
    my $NODE2 = $SPMD_NODES / 2;
    $main::MPRUN2D{'npr_yz'} = "$NODE2, 2, 2, $NODE2";
  }
  #
  # For error-growth tests
  #
  if ( defined($cam->CFG->cfg("PERGRO")) && ($cam->CFG->cfg("PERGRO") =~ /TRUE/) ) {
    $main::CAMEXP{'nestep'} = -2; 
    $main::CAMEXP{'dtime'} = 1200; 
    $main::CAMEXP{'trace_gas'}  = ".false.";
    $main::CAMEXP{'readtrace'}  = ".true.";
    $main::CAMEXP{'empty_htapes'}  = ".true.";
    $main::CAMEXP{'mfilt'}      = 145;
    if ( $adiabat =~ /false/i ) {
      $main::CAMEXP{'aqua_planet'} = ".true.";
    }
    if ( $adiabat =~ /true/i ) {
      $main::CAMEXP{'start_ymd'}    = 1231;   # Check that year/month end stuff works
      $main::CAMEXP{'start_tod'}    = 0;
    } else {
      $main::CAMEXP{'start_ymd'}    = undef;  # Use default date on file
      $main::CAMEXP{'start_tod'}    = undef;
    }
  }
  #
  # For SOM tests
  #
  if ( $cam->CFG->cfg("OCEANMODEL") eq "som" ) {
    $main::CAMEXP{'start_ymd'} = 115;   # Check that update to next SST value is ok
    $main::CAMEXP{'start_tod'} = 3600*24 - 1200*3;
    $main::CAMEXP{'inithist'}  = "\'MONTHLY\'";
  }
}

sub test_model_set_constits {
#
# Set number of advected constituents based on if trace_gas set or not
#
  my $cam = shift;
  my $trace_gas = shift;

  if ( $trace_gas =~ /true/ ) {
    $cam->CFG->setcfg( "PCNST", 7 ); # PCNST = number of advected constituents
    $cam->CFG->setcfg( "PNATS", 0 ); # PNATS = number of non-advected constituents
  } elsif ( $cam->CFG->cfg("PHYSICS") eq "ccm366" ) {
    $cam->CFG->setcfg( "PCNST", 1 );
    $cam->CFG->setcfg( "PNATS", 0 );
  } else {
    $cam->CFG->unsetcfg( "PCNST" );
    $cam->CFG->unsetcfg( "PNATS" );
  }
}
