
#	CAM.pm			Erik Kluzek
#
#	Perl5 Base Object module to handle the building and running
#	of the Community Atmospheric Model (CAM).
#
#	Other objects are used to extend the very basic functionality
#	of this object. The most important features will be carried
#	in the objects that extend this one.
#
#	Description of methods:
#
#	Important basic methods:
#
#	new ---------------- Constructor of CAM_run object.
#	process_args ------- Process the command line arguments sent to the main script.
#	usage -------------- Terminate and return the valid syntax for command line options.
#	setenv ------------- Set the given module env variable.
#	env ---------------- Get the value of the env variable in the module.
#	exists ------------- Check if the given env variable exists.
#
#	Basic methods to build/run the model. (callable from a script)
#
#	chdir -------------- Change the directory.
#	checkdir ----------- Check that directory is as expected.
#	build_support ------ Build a support program 
#				(such as makdep -- the dependency generator).
#	build_namelist ----- Build the model namelist.
#	parse_namelist ----- Parse a model namelist file.
#
#	Hidden and utility methods
#
#      CFG ---------------- Return the CAM_config object contained herein.
#	import_env --------- Import environment variables into the object.
#	config_file -------- Filepath to the configuration cache file.
#	condense_path ------ Condense a path compensating for the "/../".
#	OS ----------------- Return the Operating system being used.
#	env ---------------- Get or set a given Env variable used in the module.
#	Namelist ----------- Return the namelist object.
#	OS_list ------------ Return the list of OS's build system ported to.
#	duplicate ---------- Create a duplicate object with same data as current object.
#	interactive_build -- Set the options needed for building the model with 
#				interactive prompts.
#	print -------------- Print out contents of object.
#	clean -------------- Delete object and dependency files in BLDDIR and output files 
#				in EXEDIR.
#	arch --------------- Return the architecture type (same as preprocessor token used)
#	Platform ----------- Return a 3-letter description of Operating System 
#				(O/S) running on.
#	new_namelist ------- Construct the namelist again possibly with different
#				settings.
#	Namelist ----------- Return the namelist object that is part of this one.
#	check_gmake -------- Check to ensure you have GNUMAKE set to an something appropriate.
#	check_files_in_dir - Check that the appropriate files are in the given directory.
#	exec --------------- Execute the given system command.
#	do ----------------- Test if a given option is set (returns 0 or 1 (false or true).
#	do_build ----------- Return true if option to configure/build model is set.
#	do_log ------------- Return true if option to produce a run-log file is set.
#	do_archivelog ------ Return true if option to archive the run-log file is set.
#	do_clean ----------- Return true if option to clean old directories is set.
#	do_build_namelist -- Return true if option to build a namelist is set.
#	do_resubmit -------- Return true if option to resubmit run is set.
#	die ---------------- Do a clean termination of the script.
#                            (The primary purpose of this method is to allow
#                             it to be over-ridden or extended by other objects)
#
#	$Id: CAM.pm,v 1.17.4.11 2003/06/16 14:54:28 erik Exp $
#
use 5.004;   # Use at least this version of perl
use strict;
#use diagnostics;
use Cwd;

package CAM;
use strict;
use CAM_namelist;
use CAM_config;
use CAM_lab;
use English;

sub new {
#
# Constructor
#
  my $class = shift;
  my $self = {};

  $self->{'os'} = $OSNAME;                         # Operating system running on
  $self->{'dynamics_list'} = ["eul", "sld", "fv"]; # Valid dynamics types
  $self->{'namelist_obj'} = undef;                 # Namelist object
  $self->{'build'} = "yes";                        # If should configure and build or not
  $self->{'clean'} = "no";                         # If should clean old dirs
  $self->{'build_namelist'} = "yes";               # Build model namelist or not
  $self->{'log'} = "yes";                          # Create logfile or not
  $self->{'archivelog'} = "yes";                   # Archive logfile or not
  $self->{'resub'} = "no";                         # Resubmit simulation or not
  $self->{'logfile'} = undef;                      # Name of logfile
  $self->{'cam_lab_obj'} = CAM_lab->new;           # CAM_lab default env settings object
  #
  # Environment variables that can be set
  # (Regular data is in lower case, ENV variables are upper-case)
  #
  my @env_list = (
          "CASE",                  # Case identifier
          "CSMDATA",               # data directory
          "GNUMAKE",               # path to GNU-Make
          "RUNTYPE",               # Type of run (initial, restart, branch)
          "NCDATA_VERS",           # Version number of ncdata Initial condition files
          "LAB",                   # lab working at
          "RESUB_YEAR",            # Year to keep submitting until
          "SCRIPT_DIR"             # Directory location of script
  );
  $self->{'env_list'} = \@env_list;
  foreach my $env ( @env_list ) {
    $self->{$env} = undef;
  }
  $self->{CFG} = CAM_config->new;
  bless( $self, $class );
  $self->{config_file} = "config.".$self->Platform.".cache.xml";
  $self->import_env;
  return( $self );
}

sub arch {
#
# Return appropropriate architecture name
#
  my $self = shift;

  my $ARCH = `uname -s`;
  chomp( $ARCH );
  return( $ARCH );
}

sub condense_path {
#
# Condense the pathname compensating for "/../".
#
  my $path = shift;

  my @dirs = split "/", $path;
  my @newdirs;
  foreach my $dir (@dirs) {
    if ( $dir !~ /\.\./ ) {
        push @newdirs, $dir;
    } else {
        pop @newdirs;   # remove previous dir when current dir is..
    }
  }
  my $cpath = join '/', @newdirs;
  return( $cpath );
}

sub CFG {
#
# Return config object
#
  my $self = shift;

  return( $self->{CFG} );
}

sub check_gmake {
#
# Check that GNUMAKE is set appropriately
#
  my $self = shift;

  if ( ! defined($self->{'GNUMAKE'}) ) {
    die "ERROR:: The environment variable -- GNUMAKE is not set.\n" . 
        "        GNUMAKE must be set to to a valid version of GNU make on this system.\n";
  }
  my ($GNUMAKE) = split( / /, $self->{'GNUMAKE'} );
  my $gnumakepath = undef;
  if ( $GNUMAKE =~ /^\// ) {
     $gnumakepath = $gnumakepath;     # If using absolute path
  } else {
     my @paths = split( /:/, $ENV{PATH} );
     foreach my $path ( @paths ) {
        if ( -f "$path/$GNUMAKE" ) {
          $gnumakepath = "$path/$GNUMAKE";
          last;
        }
     }
     if ( ! defined($gnumakepath) ) {
        die "ERROR:: The environment variable -- GNUMAKE is set to ". $self->{'GNUMAKE'} . "\n" . 
            "        However, this executable name can not be found in your PATH. Please, either\n" .
            "        set GNUMAKE with the full path to $GNUMAKE or change your PATH environment \n" .
            "        variable to include the path to $GNUMAKE.\n";
     }
  }
  if ( ! -x $gnumakepath ) {
    die "ERROR:: The environment variable -- GNUMAKE is set to ". $self->{'GNUMAKE'} . "\n" . 
        "        However, this is not a valid executable program. This may be because you \n" .
        "        don't have permission to execute this program or this is not a valid GNU make\n".
        "        program to use.\n";
  }
  my $result = `$GNUMAKE -v`;
  if ( $result !~ /^GNU\s*Make\s*[version]*\s*([0-9a-z.]+)/ ) {
    die "ERROR:: The environment variable -- GNUMAKE is set to ". $self->{'GNUMAKE'} . "\n" . 
        "        However, this is not a valid version of GNU Make. Set GNUMAKE to a \n" .
        "        version of GNU Make on this system (this may require installing GNU Make)\n";
  }
}

sub check_files_in_dir {
#
# Check that the given files are in the given directory
#
  my $self = shift;
  my $dir = shift;
  my $filewildcard = shift;
  my $description = shift;

  if ( ! defined($dir) || ! defined($filewildcard) ) {
    die "ERROR:: directory or file-wildcard not sent to check_files_in_dir method\n";
  }
  if ( ! $self->CFG->exists($dir) ) {
    die "ERROR:: directory($dir) sent to check_files_in_dir method not a valid CAM_config variable\n";
  }
  if ( $self->CFG->cfg($dir) eq "builtin" ) {
    return;
  }
  my $builtin = 
        "If $dir is not needed on this platform set $dir to \"builtin\". However, \n" .
        "do this only if the Makefile is configured to build without needing $dir.\n";
  if ( ! -d $self->CFG->cfg($dir) ) {
    die "ERROR:: The environment variable $dir which should be set to a directory name is \n" . 
        "        not a valid directory! Currently $dir is set to " . $self->CFG->cfg($dir) . ".\n" .
        "        Please, set $dir to the appropriate directory name.\n" .
        "        $description\n" .
        "        $builtin\n";
  }
  my ($file) = glob( $self->CFG->cfg($dir) . "/$filewildcard"  );
  if ( ! defined($file) ) {
    die "ERROR:: The environment variable $dir which should be set to the appropriate directory\n" .
        "        does not have the files expected ($filewildcard). $description\n" .
        "        Please, set $dir to the appropriate directory name.\n"  .
        "        Currently $dir is set to " . $self->CFG->cfg($dir) . ".\n" .
        "        $builtin\n";
  }
}

sub OS {
#
# Get or Return the OS
#
  my $self = shift;
  my $value = shift;

  if ( ! defined($value) ) {
    my $OS = $self->{'os'};
    return( $OS );
  } else {
    $self->{'os'} = $value;
  }
}

sub Namelist {
#
# Set the namelist object
#
  my $self = shift;

  my $NL = $self->{'namelist_obj'};
  return( $NL );
}

sub OS_list {
#
# Return the list of OS's build system is ported to.
#
  my $self = shift;

  my @list = ( "dec_osf", "linux", "solaris", "irix", "aix" );
  return( @list );
}

sub process_args {
#
# Process the input arguments to the script
#
  my $self = shift;

  print "process_args:: Process the input arguments\n";
  while ( $ARGV[0] =~ /^.+/ ) {
    $_ = shift( @ARGV );
    if ( /^-no$/ ) {     # If one of the "-no options have been set
      $_ = shift( @ARGV );
      if ( ! defined($_) ) { $self->usage; }
      if ( /co*n*f*i*g*/ ) {         # Don't build/configure the model
        $self->do_build( "no" );
      } elsif ( /na*m*e*l*i*s*t*/ ) {  # Don't build model namelist
        $self->do_build_namelist( "no" );
      } elsif ( /lo*g*/ ) {          # Don't write out a log-file
        $self->do_log( "no" );
      } else {                       # Invalid no option
        $self->usage;
      }
    } elsif ( /^-clean/ ) {          # Clean old directories
      $self->do_clean( "yes" );
    } elsif ( /^-h/ ) {              # Help message
      $self->usage;
    } elsif ( /^-l/ ) {              # Set the Lab running at
      $_ = shift( @ARGV );
      if ( ! defined($_) ) { $self->usage; }
      $self->setenv("LAB",  $_ );     
    } elsif ( /^-/ ) {               # Invalid option
      print "Invalid option $_\n\n";
      $self->usage;
    } else {                         # Anything else must be the dynamics type
      $self->CFG->setcfg("DYNAMICS",  $_ );
    }
  }
}

sub usage {
#
# Usage statement if command arguments are done correctly
#
  my $self = shift;

  my $dyn = $self->{'dynamics_list'};
  my @dyn = @$dyn;
  my $LAB = $self->env("LAB");
  my @LAB_list = $self->{cam_lab_obj}->LAB_list;
  print <<EOF;
Usage: perl $0 dynamics [options]
Where dynamics can be set to: @dyn

Options are:

	-no value = Turn off the given option
	  -no c[onfig] = Don't configure/build the model
	  -no n[amelist] = Don't build the model namelist
	  -no l[og] = Send output to screen not logfile
	-clean = Clean any old build or executable directories
	-l = Set the lab you are running at (of @LAB_list) [$LAB]
             (Also set this by setting the Environment variable LAB)
	-h = Help (this message)

Batch operation: When submitting to a batch que you need to edit the batch
      submission information at the top of the $0 script. Most
      likely the que name and possibly the number of nodes might need to be
      changed. You need to edit either the PBS, NQS, or Loadleveler sections.

      submit to batch as either

      env SCRIPT_DIR=`pwd` qsub     $0
      env SCRIPT_DIR=`pwd` llsubmit $0

EOF
  die "Terminating";
}

sub duplicate {
#
# Create a duplicate of the current object with the same data
#
  my $self = shift;

  my $new = ref($self)->new;
  foreach my $key ( keys(%$self) ) {
    $new->{$key} = $self->{$key};
  }
  return( $new );
}

sub print {
#
# Print out contents of object
#
  my $self = shift;

  my $OS = $self->OS;
  my $LAB = $self->env("LAB");
  print "Lab running at: " . $LAB . " O/S: $OS\n";
  if ( $self->do_build ) { print "DO Build the model\n"; }
  else                   { print "Do NOT Build the model\n"; }
  if ( $self->do_build_namelist ) { print "DO Build the model namelist\n"; }
  else                            { print "Do NOT Build the model namelist\n"; }
  if ( $self->do_archivelog) { print "DO archive the logfile\n"; }
  else                       { print "Do NOT archive the logfile\n"; }
  my $RESUB_YEAR = $self->{RESUB_YEAR};
  if ( $self->do_resubmit ) { print "DO resubmit model run to year:" . $RESUB_YEAR . "\n"; }
  else                      { print "Do NOT resubmit model run\n"; }
  print "Logfile: " . $self->{'logfile'} . "\n";
}
 
sub clean {
#
# Delete the object and dependency files in the MODEL_BLDDIR and the output files
# in the MODEL_EXEDIR.
#
  my $self = shift;
  my $input = shift;

  my $pwd = cwd();
  if ( defined($input) && ($input ne "alldata") ) {
    print "Clean out the directory: $input\n";
    $self->chdir( $input );
    $self->exec( "/bin/rm -rf *.[do] *.mod esmf Depends Srcfiles *.stb *.f90" );
  } elsif ( $input eq "alldata" ) {
    my $MODEL_EXEDIR = $self->CFG->cfg( "MODEL_EXEDIR" );
    if ( -d "$MODEL_EXEDIR" ) {
      print "Clean out all data files in directory: $MODEL_EXEDIR\n";
      $self->chdir( $MODEL_EXEDIR );
      my $files = "*.nc *cam?.r* *clm?.r* lsm*";
      if ( glob($files) =~ /./ ) {
        $self->exec( "/bin/rm -f $files" );
      }
    }
  } else {
    my $MODEL_BLDDIR = $self->CFG->cfg( "MODEL_BLDDIR" );
    if ( -d "$MODEL_BLDDIR" ) {
      print "Clean out the build directory: $MODEL_BLDDIR\n";
      $self->chdir( $MODEL_BLDDIR );
      $self->exec( "/bin/rm -f *.[do] *.mod" );
    }
    my $MODEL_EXEDIR = $self->CFG->cfg( "MODEL_EXEDIR" );
    if ( -d "$MODEL_EXEDIR" ) {
      print "Clean out the executable directory: $MODEL_EXEDIR\n";
      $self->chdir( $MODEL_EXEDIR );
      my $EXENAME = $self->CFG->cfg( "EXENAME" );
      $self->exec( "/bin/rm -f compile_log.* spmdstats.* $EXENAME" );
    }
  }
  $self->chdir( $pwd );
}

sub Platform {
#
# Return a three letter string describing the platform
#
  my $self = shift;

  $_ = $self->{'os'};
  my $platform;
  if ( /dec_osf/ ) {
    $platform = "osf";
  } elsif ( /linux/ ) {
    $platform = "PC";
  } elsif ( /solaris/ ) {
    $platform = "sun";
  } elsif ( /irix/ ) {
    $platform = "sgi";
  } elsif ( /aix/ ) {
    $platform = "aix";
  } else {
    $platform = "xxx";
  }
  return( $platform );
}

sub import_env {
#
# Import the Environment variables
#
  my $self = shift;

  my $env_ref  = $self->{'env_list'};
  my @env_list = @$env_ref;
  foreach my $key ( @env_list ) {
    # Only import ENV variables if they are defined and different from current setting
    if ( defined($ENV{$key}) && ($self->{$key} ne $ENV{$key}) ) {
      $self->{$key} = $ENV{$key};
      print "Import: $key as: $ENV{$key}\n";
    }
  }
  #
  # SCRIPT_DIR and LAB
  #
  if ( defined($ENV{SCRIPT_DIR}) ) { $self->{SCRIPT_DIR} = $ENV{SCRIPT_DIR}; }
  #
  # Environment variables to set
  #
  if ( ! defined($self->env("LAB")) ) { $self->setenv("LAB",  "ncar" ); }  # Default Lab
}

sub exists{
#
# Check if the given env variable exists
#
  my $self = shift;
  my $env  = shift;

  my $class = ref($self);
  if ( ! defined($env) ) {
    die "ERROR::($class) env varible not given to exists method\n";
  }
  #
  # If env variable not used in this module
  #
  if ( ! exists($self->{$env}) ) {
    return( 0 );
  } else {
    return( 1 );
  }
}

sub env {
#
# Get the given env variable used in the module.
#
  my $self = shift;
  my $env  = shift;

  my $class = ref($self);
  if ( ! defined($env) ) {
    die "ERROR::($class) env varible not given to env method\n";
  }
  #
  # If env variable not used in this module
  #
  my $value;
  if ( ! exists($self->{$env}) ) {
    $value = "NOT-USED-HERE";
  } else {
    $value = $self->{$env};
  }
  return( $value );
}

sub setenv {
#
# Set the given env variable used in the module.
#
  my $self = shift;
  my $env  = shift;
  my $value = shift;

  my $class = ref($self);
  if ( ! defined($env) ) {
    die "ERROR::($class) env varible not given to setenv method\n";
  }
  #
  # If env variable not used in this module
  #
  if ( ! exists($self->{$env}) ) {
    die "ERROR::($class)->env env varible $env not used in this module\n";
  }
  $self->{$env} = $value;
}

sub do {
#
# Test if given option is turned on or not
# (Or set if parameter given)
#
  my $self = shift;
  my $option = shift;
  my $value = shift;

  if ( defined($value) ) {
    $self->{$option} = $value;
  } elsif ( $self->{$option} =~ /[Yy][eE]*[Ss]*/ ) {
    return( 1 );
  } else {
    return( 0 );
  }
}


sub do_build {
#
#  Test if should build/make model or not
#
  my $self = shift;
  my $value = shift;

  return( $self->do('build', $value) );
}

sub do_clean {
#
#  Test if should clean out the old directories or not.
#
  my $self = shift;
  my $value = shift;

  return( $self->do('clean', $value) );
}

sub do_build_namelist {
#
#  Test if should build model-namelist or not
#
  my $self = shift;
  my $value = shift;

  return( $self->do('build_namelist', $value) );
}

sub do_resubmit {
#
#  Test if should resubmit the model simulation or not
#
  my $self = shift;
  my $value = shift;

  return( $self->do('resub', $value) );
}

sub do_log {
#
#  Test if should create log file for run or send output to screen
#
  my $self = shift;
  my $value = shift;

  return( $self->do('log', $value) );
}

sub do_archivelog {
#
#  Test if should archive log file or not
#
  my $self = shift;
  my $value = shift;

  return( $self->do('archivelog', $value) );
}

sub chdir {
#
# Change to given directory, return error if can't
#
  my $self = shift;
  my $dir = shift;

  #
  # If directory name given is a ENV variable use it
  #
  if ( $self->exists( $dir ) ) {
    my $env = $self->env( $dir );
    $dir = $env;
  }
  #
  # If directory name given is a CAM_config variable use it
  #
  if ( $self->CFG->exists( $dir ) ) {
    my $cfg = $self->CFG->cfg( $dir );
    $dir = $cfg;
  }
  print "Change to $dir\n";
  if ( ! defined( $dir ) ) {
    die "Can not change to an empty directory";
  }
  if ( ! chdir( $dir ) ) {
    die "Could not change to directory: $dir\n";
  }
}

sub checkdir {
#
# Check that current directory is same as expected directory
#
  my $self = shift;
  my $dir = shift;
  my $die = shift;

  my $pwd = cwd( );
  if ( ! CORE::chdir( $dir ) ) {
    die "checkdir: Could not change to directory: $dir\n";
  }
  my $newpwd = cwd( );
  if ( $newpwd ne $pwd ) {
    die "$die\n\n";
  }
  if ( ! CORE::chdir( $pwd ) ) {
    die "checkdir: Could not change back to directory: $pwd\n";
  }
}

sub exec {
#
# Method to execute given command line to Bourne shell
#
  my $self = shift;
  my $command = shift;
  my $die_msg = shift;
  my $echo = shift;

  if ( ! defined($command)  ) {
    die "Command not given to exec method";
  }
  print "$command\n";
  # Use system, so that results will be echoed out, if echo option set
  if ( defined($echo) && ($echo eq "echo") ) {
    system( $command );
  # Otherwise, use backtics so that results of command will not be seen
  } else {
    `$command`;
  }
  #
  # Die on an error
  #
  if ( $? != 0 ) {
    if ( defined($die_msg) ) { print "$die_msg\n"; }
    $self->die( "ERROR:: Trouble executing: $command" );
    return;
  }
}

sub build_support {
#
# Build a support program needed for the process
#
  my $self = shift;
  my @Env_list = @_;

  my $nm = ref($self) . ":build_support";
  print "build_support:: Build a support program in directory:" . cwd() . "\n";
  if ( $self->do_clean ) {
    $self->exec( $self->{GNUMAKE} . " clean" );
  }
#
# Env or config variables commonly used in building support programs
#
  my $settings = undef;
  foreach my $env ( @Env_list ) {
    my $value = $self->env($env);
    if ( (! $self->exists($env)) && (! $self->CFG->exists($env)) ) {
      die "ERROR($nm) Env or CFG variable not associated with object\n\n";
    } elsif ( $self->exists($env) && ($value =~ /./) ) {
      $settings = "$env=$value $settings";
    } elsif ( $self->CFG->exists($env) && ($self->CFG->cfg($env) =~ /./) ) {
      $value = $self->CFG->cfg($env);
      $settings = "$env=$value $settings";
    }
  }
#-----------------------------------------------------------------------
# Build the executable
#-----------------------------------------------------------------------
  $self->exec( $self->{GNUMAKE}." $settings " );
}

sub new_namelist {
#
# Create a new namelist object
#
  my $self = shift;
  my $reset = shift;
  my $inter = shift;

  my $nm = ref($self)."\:\:build_namelist";
  my $nl = $self->{'namelist_obj'};
  if ( ! defined($nl) || (defined($reset) && ($reset == 1) ) ) {
    # Command-line option processing:
    my $RUNTYPE = $self->env( "RUNTYPE" );
    my %opts;
    $opts{runtype} = $RUNTYPE;
    $opts{'config'}      = $self->CFG->cfg("MODEL_BLDDIR") . "/" . $self->config_file;
    if ( ! defined($self->{'NAMELIST'}) ) {
      $self->{'NAMELIST'} = "atm.parm.$RUNTYPE";
    }
    $opts{'out'}         = $self->{'NAMELIST'};
    $opts{'case'}        = $self->env("CASE");
    $opts{'csmdata'}     = $self->{'CSMDATA'};
    $opts{'printlev'}    = 1;
    $opts{'test'}        = 1;
    $opts{'ncdata_vers'} = $self->env( "NCDATA_VERS" );
    if ( defined($inter) ) {
      $opts{'interactive'} = $inter;
    } else {
      $opts{'interactive'} = 0;
    }
    $self->{'namelist_obj'} = CAM_namelist->new( \%opts );
    $self->{'namelist_obj'}->set_namelists;
  }
}

sub config_file {
#
# Return or set the config-cache.xml filename
#
  my $self = shift;
  my $file = shift;
  
  if ( defined($file) ) {
    $self->{config_file} = $file;
  } else {
    $file = $self->{config_file};
    return( $file );
  }
}

sub build_namelist {
#
# Create the namelist based on the OS ENV variables and the values set
# in the CAMEXP and LSMEXP associative arrays
#
  my $self = shift;
  my $reset = shift;
  my $inter = shift;

  my $class = ref($self);
  if ( $self->do_build_namelist ) {
    print "build_namelist:: Build the model namelist\n";
    $self->new_namelist( $reset, $inter );
    $self->{'namelist_obj'}->build;
  }
}

sub parse_namelist {
#
# Parse a namelist file
#
  my $self = shift;

  print "parse_namelist:: Parse a model namelist file\n";
  $self->new_namelist;
  $self->{'namelist_obj'}->parse;
}

sub die {
#
# Terminate gracifully if found an error
#
  my $self = shift;
  my $desc = shift;

  die "$desc";
}

1   # To make use or require happy
