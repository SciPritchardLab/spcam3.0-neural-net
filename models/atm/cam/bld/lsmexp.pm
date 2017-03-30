#
#	lsmexp.pm			Erik Kluzek
#
#	Perl module to create a namelist for LSM.
#
#	Description of methods:
#
#	new ------------------ Constructor
#	set_output_values ---- Set output values based on precedence of the various input
#                              values and ensure that a valid namelist is produced.
#
#------------------------------------------------------------------------

use strict;
#use diagnostics;
use Cwd;

package lsmexp;

use atmlndnl;
@lsmexp::ISA = qw(atmlndnl  namelist);

#============================================================================

sub new {
#
# Constructor
#
  my $class = shift;
  my $optsref = shift;
  my $CAM_config = shift;

  my $interactive = $$optsref{'interactive'};
  my $file = $$optsref{'out'};
  my $printlev = $$optsref{'printlev'};

  my $default_vals = {};

  my $self = $class->SUPER::new( "lsmexp", \%main::LSMEXP, $interactive, $file,
				 "DefaultLSMEXPNamelist.xml", $default_vals,
                                 $CAM_config, $printlev );

  $self->{'printlev'} = $printlev;
  $self->{'optsref'}  = $optsref;

  $self->{DYNAMICS}   = $CAM_config->cfg("DYNAMICS");       # Dynamics to use
  $self->{PERGRO}     = $CAM_config->cfg("PERGRO");         # Error-growth option
  $self->{RESOLUTION} = $CAM_config->cfg("RESOLUTION");     # horizontal resolution

  bless( $self, $class );
  return( $self );
}

#============================================================================

sub set_output_values {

# Set the LSM namelist variables.

  my ($self, %settings) = @_;

  my $runtype = $settings{RUNTYPE};
  my $class = ref($self);
  my $nm = "$class\:\:set_default_values";

  my $NLref = $self->{'NLREF'};
  my $optsref = $self->{'optsref'};
  my $default_vals = $self->{'default_vals'};
  my $opt;

  # Get the default values from the XML file
  $self->get_default_values( %settings );

  # Root directory for default initial and boundary datasets
  my $rootdir;
  if (defined($optsref->{'csmdata'})) {
      $rootdir = $optsref->{'csmdata'};
  } elsif (defined $ENV{'CSMDATA'}) {
      $rootdir = $ENV{'CSMDATA'};
  } else {
      $rootdir = $default_vals->{'csmdata'};
  }
  my $datdir = "$rootdir/lnd/lsm1";

  # Initial conditions
  unless ( defined($NLref->{'finidat'}) ) {
      if ( defined($default_vals->{'finidat'}) ) {
	  $NLref->{'finidat'} = namelist::quote_string("$datdir/$default_vals->{'finidat'}");
      } else {
	  $NLref->{'finidat'} = "\'arbitrary initialization\'";
      }
  }
  if ( defined $NLref->{'finidat'} and $optsref->{'test'} ) {
      $self->checkinputfile('finidat');
  }

  # Directory containing datasets for creating initial conditions for
  # arbitrary resolution, i.e., when not setting finidat.
  $NLref->{'mksrfdir'} = "\'$datdir\'";

  # pergro
  unless ( defined($NLref->{'pergro'}) ) {
     if ( defined($default_vals->{'pergro'}) ) {
         $NLref->{'pergro'} = $default_vals->{'pergro'};
     }
  }

  # lsmgeo
  unless ( defined($NLref->{'lsmgeo'}) ) {
     if ( defined($default_vals->{'lsmgeo'}) ) {
         $NLref->{'lsmgeo'} = $default_vals->{'lsmgeo'};
     }
  }

}

#============================================================================

1   # to make use or require happy
