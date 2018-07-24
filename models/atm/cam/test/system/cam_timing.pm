#
#	cam_timing.pm			Erik Kluzek
#
#	Perl5 Object module to handle doing simple parsing of
#	the timing.* files output by the Community Atmospheric 
#      Model (CAM).
#
#	Methods:
#
#	$Id: cam_timing.pm,v 1.1.2.2 2003/03/12 20:38:01 hender Exp $
#
use 5.004;   # Use at least this version of perl
use strict;
#use diagnostics;
use Cwd;

package cam_timing;

sub new {
#
# Constructor
#
  my $class = shift;
  my $opts = shift;  # Options

  my $MODEL_EXEDIR = $$opts{dir}; # Required argument of model executable directory
  my $logfile = $$opts{log};      # Optional argument of output logfile name
  my $CASE = $$opts{case};        # Optional argument of caseid
  my $list = $$opts{list};        # Optional argument of list of timers to check
  my $nofail = $$opts{nofail};    # Optional argument of nofail option

  my $nm = $class . "::new";
  my $self = {};
  $self->{'MODEL_EXEDIR'} = $MODEL_EXEDIR; # Model executable directory
  if ( ! defined($MODEL_EXEDIR) ) {
    die "ERROR::($nm) directory name not sent to cam_timing constructor\n";
  }
  if ( ! -d "$MODEL_EXEDIR" ) {
    die "ERROR::($nm) directory name sent to cam_timing constructor not a directory\n";
  }
  $self->{'CASE'}         = $CASE;         # Case name
  $self->{'logfile'}      = $logfile;      # Save the name of the log-file
  $self->{'time'}         = undef;         # Save the maximum run time over each task/thread
  $self->{'total_time'}   = undef;         # Save the total run time for this case
  $self->{'sync'}         = undef;         # Synchronization time
  $self->{'beg_ts'}       = undef;         # Beginning simulation time-step
  $self->{'ts'}           = undef;         # Ending simulation time-step
  $self->{'tasks'}        = undef;         # Number of tasks
  $self->{'threads'}      = undef;         # Number of threads
  $self->{'pes'}          = undef;         # Number of processors
  if ( ! defined($list) ) {
    $self->{'timer_list'} = undef;         # Don't check any extra timers
  } else {
    $self->{'timer_list'} = $list;         # Timers to check
  }
  $self->{'nofail'}       = $nofail;       # No fail option
  bless( $self, $class );
  return( $self );
}

sub die {
#
# Method to die if nofail not set, or return error if set.
#
  my $self = shift;
  my $prompt = shift;

  if ( ! defined($self->{'nofail'}) || (! $self->{'nofail'} ) ) {
    die $prompt;
  } else {
    print "\n\nERROR:: $prompt \n\n Continuing because nofail is on\n\n";
  }
}

sub report_perf {
#
# Report a summary of the performance of this simulation
#
  my $self = shift;

  my $nm = ref($self) . "::report_perf";
  my $logfile = $self->{'logfile'};
  my $CASE = $self->{'CASE'};
  my $MODEL_EXEDIR = $self->{'MODEL_EXEDIR'};
  if ( defined($CASE) ) {
    print "Report performance of $CASE in directory $MODEL_EXEDIR\n\n";
  } else {
    print "Report performance of directory: $MODEL_EXEDIR\n\n";
  }
  if ( ! defined( $logfile) ) {
    print "Warning::($nm) logfile not set when report_perf called\n";
  } else {
    print "Information on timing: \n";
    if ( ! open( LOG, "<$logfile" ) ) { 
       $self->die( "ERROR::($nm) Can not open $logfile\n" ); 
       return( 1 ); 
    }
    my $complete = undef;
    $self->{'beg_ts'} = undef;
    $self->{'ts'}     = 0;
    while( $_ = <LOG> ) {
       #
       # Store model time-step as printed out to log file for history
       #
       if ( /Current step number:\s*([0-9]+)/ ) {
         $self->{'beg_ts'} = $1;
       }
       if ( /Number of completed timesteps:\s*([0-9]+)/ ) {
         $self->{'ts'} = $1;
       }
       #
       # Check that ran to completion
       #
       if ( /END OF MODEL RUN/ ) {
         $complete = 1;
         print "Simulation ran to completion here's the end of the log-file\n";
         print $_;
         while( $_ = <LOG> ) {
           print $_;
         }
       }
    }
    if ( ! defined( $complete ) ) {
      $self->die( "ERROR::($nm) Simulation did not run to completion\n" );
      return( 1 );
    }
  }
  #
  # Parse timing files
  #
  my $MODEL_EXEDIR = $self->{MODEL_EXEDIR};
  if ( ! -d $MODEL_EXEDIR ) {
     $self->die( "ERROR::($nm) MODEL_EXEDIR config var not set to valid directory = $MODEL_EXEDIR\n\n" );
     return( 1 );
  }
  my @timingfiles = glob( "$MODEL_EXEDIR/timing.*" );
  if ( $#timingfiles == -1 ) {
     $self->die( "ERROR::($nm) No timing files found in $MODEL_EXEDIR\n\n" );
     return( 1 );
  }
  if ( ! open( FILE, "<$timingfiles[0]" ) ) { 
     $self->die( "ERROR::($nm) Can not open $timingfiles[0]\n" ); 
     return( 1 ); 
  }
  $self->{'threads'} = 0;
  while( $_ = <FILE> ) {
    if ( /Name\s+Called/ ) {
      #print $_;
    }
    if ( /Stats for thread ([0-9]+):/ ) {
      $self->{'threads'} = $1;
    }
  }
  close( FILE );
  $self->{'threads'} += 1;
  $self->{'time'} = 0.0;
  $self->{'total_time'} = 0.0;
  $self->{'sync'} = 0.0;
  my %timer;
  if ( defined($self->{timer_list}) ) {
    my $list = $self->{timer_list};
    foreach my $timer ( @$list ) {
      $timer{$timer} = 0.0;
    }
  }
  foreach my $file ( @timingfiles ) {
    if ( ! open( FILE, "<$file" ) ) { 
      $self->die( "ERROR::($nm) Can not open $file\n" ); 
      return( 1 ); 
    }
    while( $_ = <FILE> ) {
      # Track maximum time spent in any task/thread (and total time)
      if ( /stepon\s+([0-9]+)\s+([0-9.]+)\s+([0-9.]+)\s+/ ) {
        $self->{'total_time'} += $2;
        if ( $2 > $self->{'time'} ) {
          $self->{'time'} = $2;
        }
      }
      # Keep track of max time for synchronization points
      if ( /sync[a-zA-Z0-9_]*\s+([0-9]+)\s+([0-9.]+)\s+([0-9.]+)\s+/ ) {
        if ( $2 > $self->{'sync'} ) {
          $self->{'sync'} = $2;
        }
      }
      # Track maximum time spent in list of timers to check
      if ( defined($self->{timer_list}) ) {
        my $list = $self->{timer_list};
        foreach my $timer ( @$list ) {
          if ( /$timer\s+([0-9]+)\s+([0-9.]+)\s+([0-9.]+)\s+/ ) {
            if ( $2 > $timer{$timer} ) {
              $timer{$timer} = $2;
            }
          }
        }
      }
    }
    close( FILE );
  }
  $self->{'tasks'} = $#timingfiles + 1;
  $self->{'pes'} = $self->{'tasks'}*$self->{'threads'};
  print "     Time: ". $self->{'time'} . " Total Time: " . $self->{'total_time'} . 
        " Tasks: " . $self->{'tasks'} . " Threads: " . $self->{'threads'};
  if ( defined( $self->{'beg_ts'} ) ) {
    print " start time-step: " . $self->{'beg_ts'};
  }
  if ( defined( $self->{'ts'} ) ) {
    print " end time-step: " . $self->{'ts'};
  }
  if ( $self->{'sync'} > 0.0 ) {
    print " sync time: " . $self->{'sync'};
  }
  if ( defined($self->{timer_list}) ) {
    print "\n";
    my $list = $self->{timer_list};
    foreach my $timer ( @$list ) {
      print " $timer time: " . $timer{$timer};
    }
  }
  print "\n\n";
  return( 0 );
}

sub compare_perf {
#
# Compare the performance of two simulations
#
  my $self = shift;
  my $cont = shift;  # Control case to compare performance to

  my $nm = ref($self) . "::compare_perf";
  my $time1 = $self->{'time'};
  my $time2 = $cont->{'time'};
  my $sync1 = $self->{'sync'};
  my $sync2 = $cont->{'sync'};

  if ( !defined($time1) ) {
    $self->die( "ERROR::($nm) report_perf not called before compare_perf\n" );
    return( 0.0 );
  }
  if ( !defined($time2) ) {
    $self->die( "ERROR::($nm) report_perf not called on object sent to compare_perf\n" );
    return( 0.0 );
  }
  my $case1 = $self->{CASE};
  if ( !defined($case1) ) {
    $case1 = $self->{MODEL_EXEDIR};
  }
  my $case2 = $cont->{CASE};
  if ( !defined($case2) ) {
    $case2 = $cont->{MODEL_EXEDIR};
  }
  my $rate = 0.0;
  if ( $time1 > 0.0 ) {
    $rate = 100.0*($time2-$time1)/$time2;
    if ( $rate < 0.0 ) {
      printf "%s %5.0f %s", "$case1 is slower than $case2 by ", -$rate, "%\n";
    } else {
      printf "%s %5.0f %s", "$case1 is faster than $case2 by ", $rate, "%\n";
    }
    if ( $self->{'pes'} != $cont->{'pes'} ) {
       $self->die( "ERROR::($nm) Number of processors is different between $case1 and $case2\n" );
    }
    if ( defined($self->{'beg_ts'}) && defined($cont->{'beg_ts'}) 
    && $self->{'beg_ts'} != $cont->{'beg_ts'} ) {
       $self->die( "ERROR::($nm) Starting model time-step is different between $case1 and $case2\n");
    }
    if ( defined($self->{'ts'}) && defined($cont->{'ts'}) 
    && $self->{'ts'} != $cont->{'ts'} ) {
       $self->die( "ERROR::($nm) Number of model time-step ran is different between $case1 and $case2\n" );
    }
    if ( $rate < -15.0 ) {
       print "$case1 is significantly slower than $case2\n\n";
    } elsif ( $rate > 15.0 ) {
       print "$case1 is significantly faster than $case2\n\n";
    }
  }
  if ( ($sync1 > 0.0) && ($sync2 > 0.0) ) {
    my $sync_rate = 100.0*($sync2-$sync1)/$sync2;
    if ( $rate < 0.0 ) {
      printf "%s %5.0f %s", "$case1 synchronization is shorter than $case2 by ", -$sync_rate, "%\n";
    } else {
      printf "%s %5.0f %s", "$case1 synchronization is longer than $case2 by ", $sync_rate, "%\n";
    }
  }
  return( $rate );
}


1   # To make use or require happy
