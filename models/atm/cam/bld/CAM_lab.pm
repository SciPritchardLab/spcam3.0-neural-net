#
#	CAM_lab.pm			Erik Kluzek
#
#	Object with default settings for ENV variables 
#      needed to run the scripts at different labs.
#
#	Methods:
#
#	new -------- Constructor.
#	default ---- Get default for a specific env variable.
#	LAB_list --- Get list of default labs.
#	def_list --- Get list of environment variables that have defaults.
#	initialize - Create the list of defaults.
#
#	$Id: CAM_lab.pm,v 1.15.6.15 2003/09/02 18:41:13 hender Exp $
#
use strict;
#use diagnostics;


package CAM_lab;
use lab_default;
use Env qw(LOGNAME);
#
# Constructor
#
sub new {
  my $class = shift;

  my $self = {};
  my @LAB_list = ("ncar", "ornl", "dao", "nersc", "llnl", "default");
  $self->{lab_list} = \@LAB_list;
  bless( $self, $class );
  $self->initialize;
  return( $self );
}

sub def_list {
#
# Get list of env variables that have defaults
#
  my $self = shift;

  my @keys = keys(%$self);
  my @env;
  foreach my $i ( @keys ) {
    if ( $i ne "lab_list" ) {
      push( @env, $i );
    }
  }
  return( @env );
}

sub default {
#
# Get the default value for the given variable
#
  my $self = shift;
  my $env = shift; # Name of env variable to get the default of
  my $CAM = shift; # CAM object

  my $class = ref($self);
  if ( ! defined($env) ) {
    die "ERROR:: Did not define a env variable to get a default for\n";
  }
  if ( ! exists($self->{$env}) ) {
    die "ERROR:: Default env variable $env does not exist in $class\n";
  }
  if ( ! defined($CAM) && (ref($CAM) =~ /CAM/) ) {
    die "ERROR:: Did not define a CAM object to the $class default method\n";
  }
  my $LAB = $CAM->{LAB}; 
  my $os = $CAM->{os};
  my $value = $self->{$env}->value($LAB, $os, $CAM);
}

sub LAB_list {
#
# Return the list of labs you can operate on
#
  my $self = shift;

  my $list = $self->{lab_list};
  return( @$list );
}

sub initialize {
#
# Location for case directory
#
  my $self = shift;

  my %CASE_DIR_ncar;
  my %CASE_DIR_ornl;
  my %CASE_DIR_dao;
  my %CASE_DIR_nersc;
  my %CASE_DIR_llnl;
  my %CASE_DIR_default;
  $CASE_DIR_ncar{'default'}    = "/ptmp/$LOGNAME";
  $CASE_DIR_ncar{'solaris'}    = "/net/sanitas/export/scratch/$LOGNAME";
  $CASE_DIR_ncar{'linux'}    = "/scratch/cluster/$LOGNAME";
  $CASE_DIR_default{'default'} = "\$MODEL_CFGDIR";
  $CASE_DIR_default{'dec_osf'} = "/tmp/$LOGNAME";
  $CASE_DIR_llnl{'aix'} = "/p/gf1/$LOGNAME";
  $CASE_DIR_llnl{'default'} = "/nfs/tmp0/$LOGNAME";
  $CASE_DIR_ornl{'dec_osf'} = "/usr/home/$LOGNAME";
  $CASE_DIR_ornl{'default'} = "/tmp/gpfs300a/$LOGNAME";
  $CASE_DIR_dao{'default'} = "/scratch/$LOGNAME";
  $CASE_DIR_nersc{'default'} = "/scratch/scratchdirs/$LOGNAME";
  $self->{CASE_DIR} = lab_default->new( $self->{lab_list}, (
                                                'ncar', \%CASE_DIR_ncar, 
                                                'ornl', \%CASE_DIR_ornl,
                                                'dao', \%CASE_DIR_dao,
                                                'llnl', \%CASE_DIR_llnl,
                                                'nersc', \%CASE_DIR_nersc,
                                                'default', \%CASE_DIR_default
                                              ) );
#
# Location to build
#
  my %BUILD_DIR_nersc;
  my %BUILD_DIR_ncar;
  my %BUILD_DIR_ornl;
  my %BUILD_DIR_llnl;
  my %BUILD_DIR_default;
  $BUILD_DIR_default{'aix'} = "\${SCRIPT_DIR}";
  $BUILD_DIR_default{'default'} = "\${CASE_DIR}";
  $BUILD_DIR_ncar{'default'} = "\${CASE_DIR}";
  $BUILD_DIR_ornl{'default'} = "\${CASE_DIR}";
  $BUILD_DIR_llnl{'default'} = "\${CASE_DIR}";
  $BUILD_DIR_nersc{'aix'} = "\${CASE_DIR}";
  $BUILD_DIR_nersc{'default'} = "\${SCRIPT_DIR}";
  $self->{BUILD_DIR} = lab_default->new( $self->{lab_list}, (
                                                'ncar', \%BUILD_DIR_ncar,
                                                'ornl', \%BUILD_DIR_ornl,
                                                'llnl', \%BUILD_DIR_llnl,
                                                'nersc', \%BUILD_DIR_nersc,
                                                'default', \%BUILD_DIR_default
                                              ) );
#
# Location to put log-files
#
  my %LOG_DIR_nersc;
  my %LOG_DIR_default;
  $LOG_DIR_default{'default'} = "\${MODEL_EXEDIR}";
  $LOG_DIR_nersc{'default'} = "\${MODEL_EXEDIR}";
  $self->{LOG_DIR} = lab_default->new( $self->{lab_list}, (
                                                'nersc', \%LOG_DIR_nersc,
                                                'default', \%LOG_DIR_default
                                              ) );
#
# Location to build model
#
  my %MODEL_BLDDIR_default;
  my %MODEL_BLDDIR_nersc;
  my %MODEL_BLDDIR_ncar;
  $MODEL_BLDDIR_default{'aix'} = "\${BUILD_DIR}/\${CASE}/obj";
  $MODEL_BLDDIR_default{'default'} = "\${BUILD_DIR}/\${CASE}/obj";
  $MODEL_BLDDIR_ncar{'default'} = "\${BUILD_DIR}/\${CASE}/obj";
  $MODEL_BLDDIR_nersc{'default'} = "\${BUILD_DIR}/\${CASE}/obj";
  $self->{MODEL_BLDDIR} = lab_default->new( $self->{lab_list}, (
                                                'ncar', \%MODEL_BLDDIR_ncar,
                                                'nersc', \%MODEL_BLDDIR_nersc,
                                                'default', \%MODEL_BLDDIR_default
                                              ) );
#
# Whether to run in SPMD mode or not
#
  my %SPMD_default;
  my %SPMD_ncar;
  my %SPMD_ornl;
  $SPMD_default{'default'} = "TRUE";
  $SPMD_default{'linux'} = "FALSE";
  $SPMD_default{'irix'} = "FALSE";
  $SPMD_ncar{'default'} = "TRUE";
  $SPMD_ncar{'dec_osf'} = "FALSE";
  $SPMD_ncar{'linux'} = "FALSE";
  $SPMD_ncar{'irix'} = "FALSE";
  $SPMD_ornl{'default'} = "TRUE";
  $SPMD_ornl{'dec_osf'} = "TRUE";
  $self->{SPMD} = lab_default->new( $self->{lab_list}, (
                                                'ornl', \%SPMD_ornl,
                                                'ncar', \%SPMD_ncar,
                                                'default', \%SPMD_default
                                              ) );
#
# Whether to run in SMP mode or not
#
  my %SMP_default;
  $SMP_default{'default'} = "TRUE";
  $self->{SMP} = lab_default->new( $self->{lab_list}, (
                                                'default', \%SMP_default
                                              ) );
#
# Name of the GNU-make command and the number of parallel processes to use
#
  my %GNUMAKE_dao;
  my %GNUMAKE_ncar;
  my %GNUMAKE_ornl;
  my %GNUMAKE_nersc;
  my %GNUMAKE_default;
  $GNUMAKE_ncar{'default'} = "gmake -j 2";
  $GNUMAKE_ncar{'solaris'} = "gmake ";
  $GNUMAKE_ncar{'irix'} = "gmake -j 8";
  $GNUMAKE_ncar{'dec_osf'} = "gmake -j 4";
  $GNUMAKE_ncar{'aix'} = "gmake -j 4";
  $GNUMAKE_ornl{'default'} = "gmake";
  $GNUMAKE_nersc{'default'} = "gmake -j 16";
  $GNUMAKE_dao{'default'} = "gmake";
  $GNUMAKE_dao{'aix'} = "gmake -j 8";
  $GNUMAKE_dao{'irix'} = "gmake -j 8";
  $GNUMAKE_dao{'linux'} = "gmake -j 2";
  $GNUMAKE_default{'default'} = "gmake";
  $self->{GNUMAKE} = lab_default->new( $self->{lab_list}, (
                                                 'dao', \%GNUMAKE_dao, 
                                                 'ncar', \%GNUMAKE_ncar, 
                                                 'ornl', \%GNUMAKE_ornl, 
                                                 'nersc', \%GNUMAKE_nersc, 
                                                 'default', \%GNUMAKE_default
                                            ) );
#
# Number of shared memory CPU's to run with
# (If set to 0, will unset these to use default values)
#
  my %SHMEM_CPUS_dao;
  my %SHMEM_CPUS_ncar;
  my %SHMEM_CPUS_ornl;
  my %SHMEM_CPUS_llnl;
  my %SHMEM_CPUS_default;
  $SHMEM_CPUS_ncar{'default'} = 2;
  $SHMEM_CPUS_ncar{'irix'}    = 4;
  $SHMEM_CPUS_ncar{'linux'}   = 2;
  $SHMEM_CPUS_ncar{'aix'}     = 0;
  $SHMEM_CPUS_ncar{'dec_osf'} = 4;
  $SHMEM_CPUS_dao{'default'}  = 4;
  $SHMEM_CPUS_dao{'irix'}     = 4;
  $SHMEM_CPUS_dao{'linux'}    = 2;
  $SHMEM_CPUS_dao{'aix'}      = 0;
  $SHMEM_CPUS_llnl{'default'} = 2;
  $SHMEM_CPUS_llnl{'aix'}     = 4;
  $SHMEM_CPUS_ornl{'default'} = 2;
  $SHMEM_CPUS_ornl{'dec_osf'} = 4;
  $SHMEM_CPUS_default{'default'} = 2;
  $SHMEM_CPUS_default{'aix'}     = 0;
  $self->{SHMEM_CPUS} = lab_default->new( $self->{lab_list}, (
                                                   'dao', \%SHMEM_CPUS_dao, 
                                                   'ncar', \%SHMEM_CPUS_ncar, 
                                                   'ornl', \%SHMEM_CPUS_ornl, 
                                                   'llnl', \%SHMEM_CPUS_llnl, 
                                                   'default', \%SHMEM_CPUS_default
                                               ) );
#
# Number of nodes when running with MPI
#
  my %SPMD_NODES_dao;
  my %SPMD_NODES_ncar;
  my %SPMD_NODES_ornl;
  my %SPMD_NODES_llnl;
  my %SPMD_NODES_default;
  $SPMD_NODES_dao{'default'}  = 2;
  $SPMD_NODES_dao{'aix'}      = 4;
  $SPMD_NODES_ncar{'default'} = 2;
  $SPMD_NODES_ncar{'aix'}     = 4;
  $SPMD_NODES_ncar{'dec_osf'} = 4;
  $SPMD_NODES_ornl{'default'} = 4;
  $SPMD_NODES_ornl{'dec_osf'} = 2;
  $SPMD_NODES_llnl{'default'} = 2;
  $SPMD_NODES_llnl{'aix'}     = 4;
  $SPMD_NODES_default{'default'} = 2;
  $self->{SPMD_NODES} = lab_default->new( $self->{lab_list}, (
                                                  'dao', \%SPMD_NODES_dao,
                                                  'ncar', \%SPMD_NODES_ncar,
                                                  'ornl', \%SPMD_NODES_ornl,
                                                  'llnl', \%SPMD_NODES_llnl,
                                                  'default', \%SPMD_NODES_default
                                               ) );
#
# SPMD Run command
#
  my %SPMD_RUNCMND_default;
  my %SPMD_RUNCMND_llnl;
  my %SPMD_RUNCMND_ornl;
  $SPMD_RUNCMND_default{'default'} = "mpirun -np \$SPMD_PROCS";
  $SPMD_RUNCMND_default{'solaris'} = "mpirun -machinefile machine -np \$SPMD_PROCS";
  $SPMD_RUNCMND_default{'aix'}     = "poe";
  $SPMD_RUNCMND_llnl{'default'}    = "mpirun -np \$SPMD_PROCS ";
  $SPMD_RUNCMND_llnl{'aix'}        = "poe";
  $SPMD_RUNCMND_ornl{'default'}    = "mpirun -np \$SPMD_PROCS ";
  $SPMD_RUNCMND_ornl{'dec_osf'}    = "prun -c \$SHMEM_CPUS -N \$SPMD_NODES " . 
                                     "-n \$SPMD_PROCS sh -c";
  $SPMD_RUNCMND_ornl{'aix'}        = "poe";
  my $whichprun = `which prun`;
  my $whichdmpirun = `which dmpirun`;
  if ( $whichprun =~ /^\/.*prun$/ ) {
    $SPMD_RUNCMND_llnl{'dec_osf'} = "prun -c \$SHMEM_CPUS -N \$SPMD_NODES " .
                                     "-n $\SPMD_PROCS sh -c";
  } elsif ( $whichdmpirun  =~ /^\/.*dmpirun$/ ) {
    $SPMD_RUNCMND_llnl{'dec_osf'}    = "dmpirun -np \$SPMD_PROCS ";
  }
  $self->{SPMD_RUNCMND} = lab_default->new( $self->{lab_list}, (
                                                  'llnl', \%SPMD_RUNCMND_llnl, 
                                                  'ornl', \%SPMD_RUNCMND_ornl, 
                                                  'default', \%SPMD_RUNCMND_default
                                               ) );
#
# NetCDF library location
#
  my %LIB_NETCDF_default;
  my %LIB_NETCDF_ncar;
  my %LIB_NETCDF_ornl;
  my %LIB_NETCDF_nersc;
  my %LIB_NETCDF_dao;
  $LIB_NETCDF_ncar{'solaris'} = "/contrib/lib";
  $LIB_NETCDF_ncar{'linux'}   = "/usr/local/netcdf/lib";
  $LIB_NETCDF_ncar{'irix'}    = "/usr/local/lib64/r4i4";
  $LIB_NETCDF_ncar{'default'} = "/usr/local/lib64/r4i4";
  $LIB_NETCDF_nersc{'default'}= $ENV{'NETCDF_DIR'} . "/lib";
  $LIB_NETCDF_ornl{'dec_osf'} = "/usr/local/lib64/r4i4";
  $LIB_NETCDF_ornl{'default'} = "/usr/local/lib64/r4i4";
  $LIB_NETCDF_dao{'default'} = "/usr/local/lib";
  $LIB_NETCDF_dao{'irix'}    = "/ford1/local/IRIX64/netcdf/lib";
  $LIB_NETCDF_dao{'linux'}   = "/usr/local/netcdf-3.4_Lahey/lib";
  $LIB_NETCDF_default{'default'} = "/usr/local/lib";
  $self->{LIB_NETCDF} = lab_default->new( $self->{lab_list}, (
                                                  'ncar', \%LIB_NETCDF_ncar, 
                                                  'ornl', \%LIB_NETCDF_ornl, 
                                                  'dao', \%LIB_NETCDF_dao, 
                                                  'nersc', \%LIB_NETCDF_nersc, 
                                                  'default', \%LIB_NETCDF_default
                                               ) );
#
# NetCDF includes location
#
  my %INC_NETCDF_default;
  my %INC_NETCDF_dao;
  my %INC_NETCDF_nersc;
  my %INC_NETCDF_ncar;
  $INC_NETCDF_ncar{'solaris'} = "/contrib/include";
  $INC_NETCDF_ncar{'linux'}   = "/usr/local/netcdf/include";
  $INC_NETCDF_ncar{'default'} = "/usr/local/include";
  $INC_NETCDF_nersc{'default'}= $ENV{'NETCDF_DIR'} . "/include";
  $INC_NETCDF_dao{'default'}  = "/usr/local/include";
  $INC_NETCDF_dao{'irix'}     = "/ford1/local/IRIX64/netcdf/include";
  $INC_NETCDF_dao{'linux'}    = "/usr/local/netcdf-3.4_Lahey/include";
  $INC_NETCDF_default{'default'} = "/usr/local/include";
  $self->{INC_NETCDF} = lab_default->new( $self->{lab_list}, (
                                                  'ncar', \%INC_NETCDF_ncar, 
                                                  'dao', \%INC_NETCDF_dao, 
                                                  'nersc', \%INC_NETCDF_nersc, 
                                                  'default', \%INC_NETCDF_default
                                               ) );
#
# MPI includes location
#
  my %INC_MPI_default;
  my %INC_MPI_dao;
  my %INC_MPI_ncar;
  my %INC_MPI_ornl;
  $INC_MPI_ncar{'solaris'} = "/contrib/include";
  $INC_MPI_ncar{'linux'}   = "/usr/local/mpich/include";
  $INC_MPI_ncar{'irix'}    = $ENV{'MPT_SGI'} . "/usr/include";
  $INC_MPI_ncar{'dec_osf'} = "/usr/include";
  $INC_MPI_ncar{'aix'}     = "builtin";
  $INC_MPI_ncar{'default'} = "/usr/local/include";
  $INC_MPI_ornl{'dec_osf'} = "builtin";
  $INC_MPI_ornl{'aix'}     = "builtin";
  $INC_MPI_ornl{'default'} = "/usr/local/include";
  $INC_MPI_dao{'default'}  = "/usr/local/include";
  $INC_MPI_dao{'irix'}     = $ENV{'MPT_SGI'} . "/usr/include";
  $INC_MPI_dao{'linux'}    = "/usr/local/mpich-1.2.1-ffc/include";
  $INC_MPI_dao{'aix'}      = "builtin";
  $INC_MPI_default{'default'} = "/usr/local/include";
  $INC_MPI_default{'aix'}  = "builtin";
  $self->{INC_MPI} = lab_default->new( $self->{lab_list}, (
                                                  'ncar', \%INC_MPI_ncar,
                                                  'ornl', \%INC_MPI_ornl,
                                                  'dao', \%INC_MPI_dao,
                                                  'default', \%INC_MPI_default
                                               ) );
#
# MPI library location
#
  my %LIB_MPI_default;
  my %LIB_MPI_ncar;
  my %LIB_MPI_ornl;
  my %LIB_MPI_dao;
  $LIB_MPI_ncar{'solaris'} = "/contrib/lib";
  $LIB_MPI_ncar{'linux'}   = "/usr/local/mpich/lib";
  $LIB_MPI_ncar{'dec_osf'} = "/usr/lib";
  $LIB_MPI_ncar{'irix'}    = $ENV{'MPT_SGI'} . "/usr/lib64";
  $LIB_MPI_ncar{'aix'}     = "builtin";
  $LIB_MPI_ncar{'default'} = "/usr/local/lib32/r4i4";
  $LIB_MPI_ornl{'dec_osf'} = "builtin";
  $LIB_MPI_ornl{'default'} = "builtin";
  $LIB_MPI_ornl{'aix'}     = "builtin";
  $LIB_MPI_dao{'default'}  = "/usr/local/lib";
  $LIB_MPI_dao{'irix'}     = $ENV{'MPT_SGI'} . "/usr/lib64";
  $LIB_MPI_dao{'linux'}    = "/usr/local/mpich-1.2.1-ffc/lib";
  $LIB_MPI_dao{'aix'}      = "builtin";
  $LIB_MPI_default{'default'} = "/usr/local/lib";
  $LIB_MPI_default{'aix'}  = "builtin";
  $self->{LIB_MPI} = lab_default->new( $self->{lab_list}, (
                                                  'ncar', \%LIB_MPI_ncar,
                                                  'ornl', \%LIB_MPI_ornl,
                                                  'dao', \%LIB_MPI_dao,
                                                  'default', \%LIB_MPI_default
                                               ) );

#
# Command to use to archive log files
# (Don't set if doesn't exist)
#
  my %ARCHIVE_ncar;
  my %ARCHIVE_default;
  $ARCHIVE_ncar{'default'} = "msrcp";
  $ARCHIVE_default{'default'} = " ";
  $self->{ARCHIVE} = lab_default->new( $self->{lab_list}, (
                                                'ncar', \%ARCHIVE_ncar, 
                                                'default', \%ARCHIVE_default
                                            ) );
#
# Batch submission command
#
  my %batch_dao;
  my %batch_ncar;
  my %batch_ornl;
  my %batch_default;
  $batch_dao{'default'} = " ";
  $batch_dao{'irix'} = "qsub";
  $batch_dao{'aix'} = "llsubmit";
  $batch_ncar{'default'} = "qsub";
  $batch_ncar{'dec_osf'} = " ";
  $batch_ncar{'aix'}     = "llsubmit";
  $batch_ornl{'aix'} = "llsubmit ";
  $batch_ornl{'default'} = " ";
  $batch_default{'default'} = " ";
  $batch_default{'aix'} = "llsubmit ";
  $self->{batch} = lab_default->new( $self->{lab_list}, (
                                             'dao', \%batch_dao, 
                                             'ncar', \%batch_ncar, 
                                             'ornl', \%batch_ornl,
                                             'default', \%batch_default
                                          ) );
}

1   # To make use or require happy
