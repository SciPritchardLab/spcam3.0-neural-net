#
#	CAM_run.pm			Erik Kluzek
#
#	Perl5 Object module to handle the building and running
#	of the Community Atmospheric Model (CAM).
#
#	The control is kept by setting of ENV variables as well
#	as CAMEXP and LSMEXP associative arrays which control
#	the namelists. By adding extra keyword value pairs in
#	the main script to these arrays new namelist items can
#	be added to the namelists.
#
#	This module extends the basic object "CAM" with important
#	useful functionality. Most likely methods that would need
#	to be edited would be in this module rather than the "CAM"
#	module.
#
#	Description of methods:
#
#	Important "hidden" methods:
#
#	run_time_env ------- Set the Env variables needed for run-time.
#	config_options ----- Set the options to send to the configure.csh script.
#	rest_pfile --------- Returns the name of the restart pointer file.
#				(Uses $main:CAMEXP{'rest_pfile'} namelist value if available)
#
#	Basic methods to build/run the model. (callable from a script)
#
#	setup_directories -- Setup the needed build and run directories, as well as 
#				Env variables.
#	machine_specs ------ Set the Env variables needed for this specific O/S.
#	configure ---------- Run the configure.csh to configure the model simulation.
#	make --------------- Build the model executable.
#	run ---------------- Run the model simulation.
#	resubmit ----------- Resubmit the model simulation.
#	msrcp -------------- Archive file on NCAR Mass Store using msrcp.
#
#	$Id: CAM_run.pm,v 1.15.4.15 2003/06/13 15:39:14 hender Exp $
#
use 5.004;   # Use at least this version of perl
use strict;
#use diagnostics;
use Cwd;

package CAM_run;
@CAM_run::ISA = "CAM";
#
# Use the shortened syntax for the following ENV variables
#
use Env qw(LOGNAME);
use CAM;
use CAM_namelist;
use CAM_lab;
use English;

sub new {
#
# Constructor (Add additional ENV variables to handle)
#
  my $class = shift;

  my $nm = "$class::new";
  my $self = $class->SUPER::new;
  my @more_env_list = (
         "ARCHIVE",                   # Name of command to archive files
         "VPATH",                     # Path for GMAKE for support programs
         "BUILD_DIR",                 # Build directory
         "CASE",                      # Case id
         "NAMELIST",                  # Name of output namelist file
         "CASE_DIR",                  # Directory for cases
         "INC_MPI",                   # MPI include directory
         "INC_NETCDF",                # NETCDF include directory
         "LIB_NETCDF",                # NETCDF library directory
         "LIB_MPI",                   # MPI library directory
         "LOG_DIR",                   # Directory to send log-files
         "SHMEM_CPUS",                # Number of shared-memory CPU's
         "SPMD",                      # Whether to run distributed or not
         "SMP",                       # Whether to run shared-memory processing or not
         "SPMD_CPUS_ON_NODE",         # Number of SPMD CPU's to use on one node
         "SPMD_NODES",                # Number of SPMD nodes to use
         "SPMD_RUNCMND"               # SPMD run command to use
  );
  my $env_ref  = $self->{'env_list'};
  my @env_list = @$env_ref;
  push( @env_list, @more_env_list );
  $self->{'env_list'} = \@env_list;
  foreach my $env ( @more_env_list ) {
    $self->{$env} = undef;
  }
  $self->{SPMD_PROCS} = undef;        # Number of SPMD processes (SPMD_NODES*SPMD_CPUS_ON_NODE)
  my %run_env_vars;
  $self->{'run_env_vars'} = \%run_env_vars;
  $self->{'CHANGELOG'} = "ChangeLog"; # Name of Change-file to add to top of log-files
  bless( $self, $class );
  $self->import_env;
  #
  # Put variables into CAM_config if defined
  #
  my @list = ("SPMD", "SMP", "INC_MPI", "INC_NETCDF", "LIB_MPI", "LIB_NETCDF" );
  foreach my $key ( @list ) {
    if ( defined($self->{$key}) && ($self->{$key} =~ /./) ) {
      $self->CFG->setcfg( $key, $self->{$key} );
    } else {
      $self->CFG->setcfg("$key",  $self->{cam_lab_obj}->default("$key", $self) ); 
    }
    # Delete from this object, so can only deal with it in the CAM_config object
    delete( $self->{$key} );
  }
  #
  # Set GNUMAKE
  #
  my $GNUMAKE = $self->env("GNUMAKE");
  if ( ! defined($GNUMAKE) ) { 
    $self->setenv("GNUMAKE",  $self->{cam_lab_obj}->default("GNUMAKE", $self) ); 
  }
  $self->check_gmake;
  #
  # Make sure library locations are set
  #
  $self->check_files_in_dir( "INC_NETCDF", "netcdf.inc",  
        "INC_NETCDF is the include directory for NetCDF library" );
  $self->check_files_in_dir( "LIB_NETCDF", "libnetcdf.a", 
        "LIB_NETCDF is the location for NetCDF library" );
  $self->check_files_in_dir( "LIB_MPI",    "libmpi*",   
        "LIB_MPI is the location for MPI library" );
  $self->check_files_in_dir( "INC_MPI",    "mpif.h",      
        "INC_MPI is the include directory for MPI library" );
  return( $self );
}

sub run_time_env {
#
# Run-time Platform dependencies...
#
  my $self = shift;

  my $os = $self->{'os'};
  my $LAB = $self->{LAB};
  my $run_ref = $self->{'run_env_vars'};

  if ( ! defined($self->{SPMD_NODES}) ) { 
     $self->{SPMD_NODES} = $self->{cam_lab_obj}->default("SPMD_NODES", $self);
  }
  # If SPMD_NODES set to zero, change to undefined so that defaults will be used
  if ( $self->{SPMD_NODES} == 0 ) { $self->{SPMD_NODES} = undef; }
  SWITCH: {
#---------------------------------------------------------------------------------------
# Solaris
#---------------------------------------------------------------------------------------
    ($os =~ /solaris/) && do {
                  #
                  # Prepare hostfilename of valid hosts to use
                  # Must agree with filename used in SPMD_RUNCMND (check CAM_lab.pm)
                  #
                  my $filename = "machine";
                  open( MACH, ">$filename" ) || die "ERROR:: Can not open file: $filename\n";
                  print MACH `hostname`;
                  close( MACH );
                  system( "ulimit -s unlimited" );
               last SWITCH;
             };
    ($os =~ /irix/) && do {
#---------------------------------------------------------------------------------------
# SGI
#---------------------------------------------------------------------------------------
# OMP_DYNAMIC       : False => Don't give up processors even if machine is
#                     busy
# _DSM_PLACEMENT    : Ensure efficient memory allocation 
# _DSM_WAIT         : Assume dedicated access to processors (SPIN)
# _DSM_VERBOSE      : Print diagnostic info about system-level software
# MPC_GANG          : Gang scheduling.  Does not work well so turn off
# MP_SLAVE_STACKSIZE: Stack size (bytes) to be used by slave processes
# TRAP_FPE          : Abort on overflow or divide by zero
#
#---------------------------------------------------------------------------------------
                  my %env = (OMP_DYNAMIC=>"FALSE", _DSM_PLACEMENT=>"ROUND_ROBIN", 
                             _DSM_WAIT=>"SPIN",    MPC_GANG=>"OFF", 
                             MP_SLAVE_STACKSIZE=>"40000000", 
                      TRAP_FPE=>"UNDERFL=FLUSH_ZERO; OVERFL=ABORT,TRACE; DIVZERO=ABORT,TRACE" );
                  foreach my $key ( keys(%env) ) {
                    if ( $ENV{$key} !~ /.+/ ) { 
                       $$run_ref{$key} = $env{$key};
                    } else {
                       $$run_ref{$key} = $ENV{$key};
                    }
                  }
                  system( "ulimit -s unlimited" );
               last SWITCH;
             };
    ($os =~ /linux/) && do {
#---------------------------------------------------------------------------------------
# Linux
#---------------------------------------------------------------------------------------
# MPSTKZ  sets slave stack size
#---------------------------------------------------------------------------------------
                  my $key = "MPSTKZ";
                  if ( $ENV{$key} !~ /.+/ ) { 
                     $$run_ref{$key} = "128M";
                  } else {
                     $$run_ref{$key} = $ENV{$key};
                  }
               last SWITCH;
             };
    ($os =~ /aix/) && do {
#---------------------------------------------------------------------------------------
# IBM
#---------------------------------------------------------------------------------------
# XLSMPOPTS             the largest amount of space (bytes) that a thread's 
#                       stack will need 
# MP_EUILIB us          message passing subsystem implementation
# MP_NODES              number of physical nodes to run the parallel tasks
#                       (a task refers specifically to an MPI process)
# MP_TASKS_PER_NODE     number of tasks to run on each of the physical nodes 
# MP_RMPOOL             specifies number of a LoadLeveler pool     
#---------------------------------------------------------------------------------------
                  my %env = ( XLSMPOPTS=>"stack=40000000", MP_EUILIB=>"us",
                              MP_RMPOOL=>"1" );
                  foreach my $key ( keys(%env) ) {
                    if ( $ENV{$key} !~ /.+/ ) { 
                       $$run_ref{$key} = $env{$key};
                    } else {
                       $$run_ref{$key} = $ENV{$key};
                    }
                  }
                  #
                  # If SPMD_NODES not set (either value or ENV var, but MP_NODES is set)
                  # use the MP_NODES value
                  #
                  if ( (! defined($self->{SPMD_NODES})) && ($ENV{'SPMD_NODES'} !~ /./) 
                  && ($ENV{'MP_NODES'} =~ /./) ) { 
                    $self->{SPMD_NODES} = $ENV{'MP_NODES'};
                  }
                  $$run_ref{MP_NODES} = $self->{SPMD_NODES};
                  #
                  # If SPMD_CPUS_ON_NODE not set, use a default of 1 or
                  # use the MP_TASKS_PER_NODE value
                  #
                  if ( (! defined($self->{SPMD_CPUS_ON_NODE})) && ($ENV{'SPMD_CPUS_ON_NODE'} !~ /./) 
                  && ($ENV{'MP_TASKS_PER_NODE'} =~ /./) ) { 
                    $self->{SPMD_CPUS_ON_NODE} = $ENV{'MP_TASKS_PER_NODE'};
                  }
                  if ( ! defined($self->{SPMD_CPUS_ON_NODE}) ) { 
                    $self->{SPMD_CPUS_ON_NODE} = 1; 
                  }
                  $$run_ref{MP_TASKS_PER_NODE} = $self->{SPMD_CPUS_ON_NODE};
               last SWITCH;
             };
    ($os =~ /dec_osf/) && do {
#---------------------------------------------------------------------------------------
# Compaq
#---------------------------------------------------------------------------------------
# MP_STACK_SIZE  slave stack size
#---------------------------------------------------------------------------------------
                  my $key = "MP_STACK_SIZE";
                  if ( $ENV{$key} !~ /.+/ ) { 
                     $$run_ref{$key} = "17000000";
                  } else {
                     $$run_ref{$key} = $ENV{$key};
                  }
                  system( "ulimit -s unlimited" );
               last SWITCH;
             };
#---------------------------------------------------------------------------------------
# default
#---------------------------------------------------------------------------------------
  }
  #
  # If SHMEM_CPUS not set (either value or ENV var, but OMP_NUM_THREADS is set)
  # use the OMP_NUM_THREADS value
  #
  if ( (! defined($self->{SHMEM_CPUS})) && ($ENV{'SHMEM_CPUS'} !~ /./) 
  && ($ENV{'OMP_NUM_THREADS'} =~ /./) ) { 
    $self->{SHMEM_CPUS} = $ENV{'OMP_NUM_THREADS'};
  }
  if ( ! defined($self->{SHMEM_CPUS}) ) { 
    $self->{SHMEM_CPUS} = $self->{cam_lab_obj}->default("SHMEM_CPUS", $self);
  }
  #
  # SPMD cpus on node
  #
  if ( ! defined($self->{SPMD_CPUS_ON_NODE}) ) { 
    $self->{SPMD_CPUS_ON_NODE} = 1; 
  }
  # SPMD processes
  $self->{SPMD_PROCS} = $self->{SPMD_NODES} * $self->{SPMD_CPUS_ON_NODE};

  # If SHMEM_CPUS set to zero, change to undefined so that defaults will be used
  if ( $self->{SHMEM_CPUS} == 0 ) { $self->{SHMEM_CPUS} = undef; }
  if ( defined($self->{SHMEM_CPUS}) ) { 
    $$run_ref{'OMP_NUM_THREADS'} = $self->{SHMEM_CPUS}; 
  }
  if ( ! defined($self->{SPMD_RUNCMND}) ) { 
     $self->{SPMD_RUNCMND} = $self->{cam_lab_obj}->default("SPMD_RUNCMND", $self);
  }
}

sub config_options {
#
# Return the options to pass to configure, based on what
# ENV variables have been set.
#
  my $self = shift;

  my $nm = ref($self) . "::config_options ";
  # Set SPMD to default setting
  my $LAB = $self->{LAB};
  print "config_options:: Figure out the options to pass to configure\n";
  my %opt = ( PCNST=>'-nadv', PNATS=>'-nnadv', RESOLUTION=>'-res',
                 PLEV=>'-nlev', DYNAMICS=>'-dyn', MODEL_BLDDIR=>'-cam_bld',
                 CAMROOT=>'-cam_root', MODEL_EXEDIR=>'-cam_exedir',
                 MODEL_CFGDIR=>'-cam_cfg', PHYSICS=>'-phys',
                 MODEL_MODDIR=>'-usr_src',
                 EXENAME=>'-cam_exe', OCEANMODEL=>'-ocn', ICEMODEL=>'-sice'
             );
  my $options = "-cache " . $self->config_file;
  foreach my $key ( keys(%opt) ) {
     if ( defined($self->CFG->cfg($key)) && ($self->CFG->cfg($key) =~ /.+/) ) {
        $options = $options . " $opt{$key} " . $self->CFG->cfg($key);
     }
  }
  # Resolution options
  %opt = ( PLAT=>'-nlat', PLON=>'-nlon', PTRM=>'-trm', PTRN=>'-trn', PTRK=>'-trk',
         );
  if ( $self->CFG->cfg("RESOLUTION") eq "custom" ) {
    foreach my $key ( keys(%opt) ) {
       if ( defined($self->CFG->cfg($key)) && ($self->CFG->cfg($key) =~ /.+/) ) {
          $options = $options . " $opt{$key} " . $self->CFG->cfg($key);
       } else {
         my @keys = keys(%opt);
         die "ERROR($nm): Resolution = \`custom\`, but needed settings aren't provided: @keys\n";
       }
    }
  }
  # Library directories
  my %eopts = ( INC_MPI=>'-mpi_inc', INC_NETCDF=>'-nc_inc', 
                LIB_MPI=>'-mpi_lib', LIB_NETCDF=>'-nc_lib' );
  foreach my $key ( keys(%eopts) ) {
     my $value = $self->CFG->cfg($key);
     if ( defined($value) && ($value =~ /./) && ($value !~ /builtin/) ) {
        $options = $options . " $eopts{$key} " . $value;
     }
  }
  if ( defined($self->CFG->cfg("SPMD")) ) {
    if ( ($self->CFG->cfg("SPMD") =~ /TRUE/) ) {
      $options = $options . " -spmd ";
     } elsif ( ($self->CFG->cfg("SPMD") =~ /FALSE/) ) {
       $options = $options . " -nospmd ";
     } else {
       die "ERROR($nm): Invalid value for SPMD = " . $self->CFG->cfg("SPMD") . "\n";
     }
  }
  if ( defined($self->CFG->cfg("SMP")) ) { 
    if ( ($self->CFG->cfg("SMP") =~ /TRUE/) ) {
      $options = $options . " -smp ";
    } elsif ( ($self->CFG->cfg("SMP") =~ /FALSE/) ) {
      $options = $options . " -nosmp ";
    } else {
       die "ERROR($nm): Invalid value for SMP = " . $self->CFG->cfg("SMP") . "\n";
    }
  }
  if ( defined($self->CFG->cfg("TWOD_YZ")) ) { 
    if ( ($self->CFG->cfg("TWOD_YZ") =~ /TRUE/) ) {
      $options = $options . " -twod_yz ";
    } elsif ( ($self->CFG->cfg("TWOD_YZ") =~ /FALSE/) ) {
    } else {
       die "ERROR($nm): Invalid value for TWOD_YZ = " . $self->CFG->cfg("TWOD_YZ") . "\n";
    }
  }
  if ( defined($self->CFG->cfg("PERGRO")) && ($self->CFG->cfg("PERGRO") =~ /TRUE/) ) {
    $options = $options . " -pergro";
  }
  if ( defined($self->CFG->cfg("DEBUG")) && ($self->CFG->cfg("DEBUG") =~ /TRUE/) ) {
    $options = $options . " -debug";
  }
  return( $options );
}

sub configure {
#
# Run the configure script
#
  my $self = shift;

  #
  # Remove Makefile if it exists
  #
  my $MODEL_BLDDIR = $self->CFG->cfg("MODEL_BLDDIR");
  if ( -l "$MODEL_BLDDIR/Makefile" ) {
    system( "/bin/rm $MODEL_BLDDIR/Makefile" );
  }
  if ( $MODEL_BLDDIR ne $self->CFG->cfg("MODEL_CFGDIR")) { 
    if ( -f "$MODEL_BLDDIR/Makefile" ) {
      system( "/bin/rm $MODEL_BLDDIR/Makefile" );
    }
  }
  #
  # Save INC_MPI and LIB_MPI in case configure changes them
  #
  my $INC_MPI = $self->CFG->cfg("INC_MPI");
  my $LIB_MPI = $self->CFG->cfg("LIB_MPI");
  # Figure out options to pass to configure based on ENV vars
  my $options = $self->config_options;
  my $config = $self->CFG->cfg("MODEL_CFGDIR")."/configure";
  if ( ! -f "$config" ) {
    die "Error:: Configure file: $config not available -- MODEL_CFGDIR not correct?";
  }
  print "$config $options\n";
  #
  # Move existing preprocessor files over to temporary names
  #
  my @files = ("misc.h", "params.h", "preproc.h");
  foreach my $file ( @files ) {
    if ( -f $file ) {
      system( "/bin/mv $file $file.tmp" );
    }
  }
  system( "$config $options" );    # Configure model 
  if ( $? != 0 ) { die "Error in configure"; }
  my $MODEL_BLDDIR = $self->CFG->cfg("MODEL_BLDDIR");
  my $config_file = "$MODEL_BLDDIR/" . $self->config_file;
  $self->CFG->read_config_cache( $config_file );
  $self->CFG->setcfg("INC_MPI", $INC_MPI);  # Reset INC_MPI back
  $self->CFG->setcfg("LIB_MPI", $LIB_MPI);  # Reset LIB_MPI back
  #
  # Move temporary preprocessor files back if same as new files created
  #
  foreach my $file ( @files ) {
    if ( -f "$file.tmp" ) {
      system( "cmp $file $file.tmp" );
      if ( $? == 0 ) {
        system( "/bin/mv -f $file.tmp $file" );
      }
    }
  }
  my $SCRIPT_DIR = $self->env( "SCRIPT_DIR" );
  system( "/bin/cp $config_file $SCRIPT_DIR" );
}

sub make {
#
# Build the model using GNU make
#
  my $self = shift;


#
# Keep track of compile log, add date of compile to log file
#
  print "make:: Make the model executable\n";
  my $log = $self->{LOG_DIR}."/compile_log.atm";
  print "Compiling CAM ... see $log for log\n\n";

  open( LOG, ">$log" ) || die "Can not add output to the compile log file: $log\n";
  print LOG << "EOF";
-------------------------------------
`date`
-------------------------------------
EOF
  close( LOG );
#-----------------------------------------------------------------------
# Build the executable
#-----------------------------------------------------------------------
  $self->exec( $self->{GNUMAKE}." >> $log 2>&1", "Compile failed look at $log for failure" );
  print "Compile successful!\n\n";
}

sub setup_directories {
#
# Setup the directories
#
  my $self = shift;

  print "setup_directories:: Setup the directories for building and running the model:\n";
  print "case name and title: ".$self->env("CASE")."\n";
  my $OS = $self->OS;
  my $LAB = $self->{'LAB'};
  # If irt set to zero, don't archive mass store files
  if ( defined($main::CAMEXP{'mss_irt'}) && ($main::CAMEXP{'mss_irt'} == 0) ) {
    $self->do_archivelog( "no" );
  }
  #
  # Case directory (directory where output from different cases are stored)
  #
  if ( ! defined($self->{CASE_DIR}) ) {
    $self->{CASE_DIR} = $self->{cam_lab_obj}->default("CASE_DIR", $self);
  }
  print "dir: ".$self->{CASE_DIR}." lab: $LAB, OS: $OS\n";
  if ( ! -d $self->{CASE_DIR} ) { 
     mkdir( $self->{CASE_DIR}, 0755 ) || die "Can not mkdir ".$self->{CASE_DIR}. "\n"; 
  }
  print "Case directory: ".$self->{CASE_DIR}."\n";
  #
  # Root directory
  #
  my $ROOT = $self->CFG->cfg("CAMROOT");
  if ( ! defined($ROOT) ) {
    $ROOT = cwd( ) . "/../../../..";
    $self->CFG->setcfg( "CAMROOT", $ROOT );
  }
  if ( ! -d $ROOT ) { 
    die "Error:: Root directory: $ROOT does not exist";
  }
  print "Root directory: $ROOT\n";
  #
  # Config directory
  #
  my $MODEL_CFGDIR = $self->CFG->cfg("MODEL_CFGDIR");
  if ( ! defined($MODEL_CFGDIR) ) {
     $self->CFG->setcfg( "MODEL_CFGDIR",  "$ROOT/models/atm/cam/bld" );
     $MODEL_CFGDIR = $self->CFG->cfg("MODEL_CFGDIR");
  }
  if ( ! -d $MODEL_CFGDIR ) { 
    die "Error:: Configuration file directory: $MODEL_CFGDIR does not exist";
  }
  print "Config directory: $MODEL_CFGDIR\n";
  #
  # Execution directory (Need CASE_DIR and CASE set to get it)
  #
  # Case name
  if ( ! defined($self->env("CASE")) ) {
    die "Error:: Env variable CASE not defined";
  }
  # Now set EXEDIR
  if ( ! defined($self->CFG->cfg("MODEL_EXEDIR")) ) {
    $self->CFG->setcfg("MODEL_EXEDIR", $self->{CASE_DIR}."/".$self->env("CASE") );
  }
  print "Exec directory: ".$self->CFG->cfg("MODEL_EXEDIR")."\n";
  #
  # Build directory
  #
  if ( ! defined($self->{BUILD_DIR}) ) {
    $self->{BUILD_DIR} = $self->{cam_lab_obj}->default("BUILD_DIR", $self);
  }
  if ( ! defined($self->CFG->cfg("MODEL_BLDDIR")) ) {
    $self->CFG->setcfg("MODEL_BLDDIR", $self->{cam_lab_obj}->default("MODEL_BLDDIR", $self) );
  }
  #
  # Directory where log-files go
  #
  $self->{LOG_DIR} = $self->{cam_lab_obj}->default("LOG_DIR", $self);
  print "dir: ".$self->{LOG_DIR}." lab: $LAB, OS: $OS\n";
  if ( ! -d $self->{LOG_DIR} ) { 
     mkdir( $self->{LOG_DIR}, 0755 ) || die "Can not mkdir ".$self->{LOG_DIR}."\n"; 
  }
  print "Log-file directory: ".$self->{LOG_DIR}."\n";
  #
  # Clean blddir and exedir if requested
  #
  if ( $self->do_clean ) {
    $self->clean;
  }
  #
  # Create exedir and blddir
  #
  if ( ! -d $self->CFG->cfg("MODEL_EXEDIR") ) { 
    mkdir( $self->CFG->cfg("MODEL_EXEDIR"), 0755 ) || die "Can not mkdir ".$self->CFG->cfg("MODEL_EXEDIR")."\n"; 
  }
  if ( ! -d $self->CFG->cfg("MODEL_BLDDIR") ) { 
    mkdir( $self->CFG->cfg("MODEL_BLDDIR"), 0755 ) || die "Can not mkdir ".$self->CFG->cfg("MODEL_BLDDIR")."\n"; 
  }
  print "Build directory: ".$self->CFG->cfg("MODEL_BLDDIR")."\n";
}

sub run {
#
# Run the model
#
  my $self = shift;
  my $desc = shift;
  my $logfile = shift;

  print "run:: Run the model\n";
  if ( defined($desc) ) { 
    print "$desc\n";
    print STDERR "$desc\n";
  } 
  #
  # Check that in $MODEL_EXEDIR
  #
  my $pwd = cwd;
  $self->checkdir( $self->CFG->cfg("MODEL_EXEDIR"), 
     "Error: Not calling run from within \$MODEL_EXEDIR(".$self->CFG->cfg("MODEL_EXEDIR").")" );
  #
  # Check that Namelist exists
  #
  my $NAMELIST = $self->env("NAMELIST");
  if ( ! -f $NAMELIST ) { die "Error: Namelist $NAMELIST does not exist"; }
  if ( -z $NAMELIST ) { die "Error: Namelist $NAMELIST is empty!"; }

  my $exec = "./".$self->CFG->cfg("EXENAME");
  if ( ! -f $exec ) {
    die "Error: Executable ($exec) does not exist!";
  }
  #
  # Check that SPMD_RUNCMND set if SPMD=TRUE
  #
  if ( $self->CFG->cfg("SPMD") eq "TRUE" && ! defined($self->env("SPMD_RUNCMND")) ) {
    die "Error: run_time_env method not called before run method";
  }
  #
  my $OS = $self->OS;
  my $log; 
  #
  # Prepare log file, put namelist and ChangeLog at top of log
  #
  if ( $self->do_log ) {
    if ( ! defined($logfile) ) {
      my $logid = `/bin/date +%y%m%d-%H%M%S`; chomp( $logid );
      $logfile = $self->{LOG_DIR}."/cam2.".$self->env("CASE").".log.$logid";
    }
    $self->{'LOGFILE'} = $logfile;
    $log = ">> $logfile 2>&1";
    if ( -f $logfile ) { unlink( $logfile ); }
    open( LOG, ">$logfile" ) || die "Can not open logfile $logfile";
    open( NAMELIST , "<$NAMELIST" ) || die "Can not open namelist $NAMELIST";
    print LOG "Namelist:\n";
    while( defined($_ = <NAMELIST>) ) { print LOG $_; }
    close( NAMELIST );
    my $change = $self->{CHANGELOG};
    if ( defined($change) ) {
      $change = $self->CFG->cfg("CAMROOT")."/models/atm/cam/doc/$change";
      $change = &CAM::condense_path( $change );
      print "Change: $change\n";
      if ( -f $change ) {
        open( CHANGELOG , "<$change" ) || die "Can not open ChangeLog $change";
        print LOG "ChangeLog:\n";
        <CHANGELOG>;
        while( defined($_ = <CHANGELOG>) && (!/====/) ) { print LOG $_; }
        close( CHANGELOG );
      }
    }
    close(LOG);
  #
  # If no-log file option used
  #
  } else {
    $log = " ";
    $logfile = " ";
  }
  print "Running CAM, ".$self->env("RUNTYPE")."... directory " . cwd() . "\n\n";
  #
  # Prepare environment variable settings that must be used
  #
  my $run_ref = $self->{'run_env_vars'};
  my $environ = "";
  foreach my $key ( keys(%$run_ref) ) {
    $environ = "$environ $key=\"" . $$run_ref{$key} . "\"";
  }
  if ( $environ ne "" ) {
    $environ = "env $environ";
  }
  #
  # Actually run the model
  #
  my $die_msg = "Error running the model";
  if ( $self->do_log ) { $die_msg = "$die_msg - see $logfile for log \n\n"; }
  if ( $self->CFG->cfg("SPMD") ne "TRUE" ) {
    $self->exec( "$environ time $exec < $NAMELIST $log", $die_msg, "echo" );
  } elsif ( $self->{SPMD_RUNCMND} =~ /sh -c/ ) {
    $self->exec( "$environ time ".$self->{SPMD_RUNCMND}." \'$exec < $NAMELIST $log\'", 
                 $die_msg, "echo" );
  } else {
    $self->exec( "$environ time ".$self->{SPMD_RUNCMND}." $exec < $NAMELIST $log", 
                 $die_msg, "echo" );
  }
  print "\n\nCAM Finished";
  if ( $self->do_log ) { print " - see $logfile for log \n\n"; }
  print "\n\n";
  #
  # Save Namelist, logfile and config-file to Mass Store system (if available)
  #
  my $LAB = $self->{LAB};
  if ( ! defined($self->{ARCHIVE}) ) {
    $self->{ARCHIVE} = $self->{cam_lab_obj}->default("ARCHIVE",$self);
  }
  if ( $self->do_archivelog && defined($self->{ARCHIVE}) && $self->do_log ) {
    #
    # Parse namelist to get MSS info
    #
    if ( ! defined($self->Namelist) ) {
      $self->new_namelist( 1 );
    }
    $self->Namelist->parse( $NAMELIST );
    $self->msrcp( "logs", $NAMELIST );
    $self->msrcp( "logs", $logfile );
    my $config_file = $self->config_file;
    if ( -f $config_file ) {
      $self->msrcp( "logs", $config_file );
    }
  }
}

sub msrcp {
#
# Executate msrcp to store a file for this case.
#
  my $self = shift;
  my $subdir = shift;         # Subdirectory (relative to ARCHIVE_DIR on namelist) to store file
  my $local_filepath = shift; # Full local disk filepath for file to store
  my $get = shift;            # Optional argument, if want to get rather than put

  my $nm = ref($self) . "::msrcp";
  #
  # Figure out full path to msrcp and make sure it exists
  #
  my $cmd = $self->{cam_lab_obj}->default("ARCHIVE",$self);
  if ( $cmd !~ /msrcp/ ) {
    return;
  }
  my $archive = `/bin/which $cmd`;   # Test if command actually exists
  chomp( $archive );
  if ( $? != 0 ) {
    print "WARNING: Archive ($cmd) does not exist in PATH\n\n";
    return;
  }
  #
  # Split out filename to the pathname and filename
  #
  my $file;
  ($file = $local_filepath) =~ s!(.*)/!!;
  my $local_dir = $1;
  if ( ! defined($local_dir) ) {
    $local_dir = ".";
  }
  #
  # Check namelist hash for retention period and archive_dir
  #
  my $reten = 365;
  if ( defined($main::CAMEXP{mss_irt}) ) {
     $reten = $main::CAMEXP{mss_irt};
  }
  my $MSSNAME = $ENV{LOGNAME};
  $MSSNAME =~ tr/[a-z]/[A-Z]/;
  my $dir = "/$MSSNAME/csm/".$self->env("CASE")."/atm";
  if ( defined($main::CAMEXP{archive_dir}) ) {
    $dir = $main::CAMEXP{archive_dir};
  }
  #
  # Finally save the file to MSS
  #
  my $cmd;
  if ( ! defined($get) || (! $get) ) {
    print "Save $file to $dir/$subdir\n";
    $cmd = "$archive -pe $reten $local_filepath mss:$dir/$subdir/$local_filepath";
  } else {
    print "Get $file from $dir/$subdir\n";
    $cmd = "$archive mss:$dir/$subdir/$file $local_dir";
  }
  print "$cmd\n";
  `$cmd`;
  if ( $? != 0 ) {
     print "ERROR::($nm):: on msrcp : $cmd";
  }
}

sub rest_pfile {
#
# Get the name for the restart pointer file.
#
  my $self = shift;

  my $file;
  my $nm = ref($self) . ":rest_pfile";
  if ( defined($main::CAMEXP{'rest_pfile'}) ) {
    $file = $main::CAMEXP{'rest_pfile'};
  } else {
    my $CASE = $self->env("CASE");
    if ( ! defined($CASE) ) {
      die "ERROR::($nm) Variable CASE not defined, can not get restart pointer\n";
    }
    $file = "$ENV{'HOME'}/cam2.$CASE.rpointer";
  }
}

sub resubmit {
#
# Resubmit the model simulation
#
  my $self = shift;

  # If resubmit option not set just return
  if ( ! $self->do_resubmit ) { 
    return;
  }
  # If not a restart run-type, print error and return
  if ( $self->env("RUNTYPE") ne "restart" ) {
    print "Trying to resubmit a simulation that is not a restart (set RUNTYPE to restart)\n";
    return;
  }
  # Read restart pointer file, keep submitting until hits RESUB_YEAR
  print "resubmit:: Resubmit this model simulation until it reaches year: ".$self->{RESUB_YEAR}."\n";
  my $file = $self->rest_pfile;
  if ( ! open( RPOINT, "<$file" ) ) {
     $self->die( "Could not open restart pointer file: $file\n" );
     return;
  }
  $_ = <RPOINT>;
  if ( /cam2.r.(\d+)-/ ) {
     my $year = $1;
     if ( $year < $self->{RESUB_YEAR} ) {
       my $LAB = $self->{LAB};
       my $OS = $self->OS;
       my $batch = $self->{cam_lab_obj}->default("batch",$self);
       if ( $batch !~ /^ $/ ) {
         # Check that in correct directory
         my $pwd = cwd( );
         my $SCRIPT_DIR = $self->{SCRIPT_DIR};
         if ( $pwd !~ /$SCRIPT_DIR/ ) { 
           print "Warning: resubmit: Need to call resubmit from $SCRIPT_DIR, currently in $pwd\n"; 
         }
         # IBM retains a ll path when scripts have been submitted to batch
         my $ScriptName;
         ($ScriptName = $PROGRAM_NAME) =~ s!(.*)/!!;         # name of program
         my $ScriptDir = $1;                      # name of directory where script lives
         $self->exec( "$batch $ScriptName", "error in batch submission", "echo" );
       } else {
         print "I do not know what the batch submission command is for: " . 
               "$LAB $OS\n";
       }
     } else {
       print "resubmit:: Beyond resubmission year:\n";
     }
   } else {
     print "resubmit:: Can not interpret the restart pointer file format\n";
   }
}

1   # To make use or require happy
