#!/usr/bin/perl
#LLNL
#/usr/local/bin/perl
package CCM_Env;
use Cwd;
use File::Basename;
use Shell qw(uname);
use Carp;

my %bool = (0 => FALSE, 1 => TRUE);
my %fields = (
    Dyn     => 1         ,
    Phy     => 1         ,
    Ver     => "fvcam"  ,
    Case    => "b26"     ,
    Obj     => "obj"     ,
    USER_FC => undef     ,
    Decomp  => twod      ,
    MPI     => 1         ,
    OMP     => 1         ,
    N_MPI   => 4         ,
    N_TOTAL => undef     ,
    NPR_YZ  => "1,1,1,1" ,
    N_SMP   => 4         ,
    EXE     => "atm"     ,
    Home    => cwd()     ,
#   FVData  => "/share/fvgcm/CAM/inputdata/atm/cam2" ,
#   LMData  => "/share/fvgcm/CAM/inputdata/lnd/clm2" ,
#   Scratch => '$SCRATCH1',
#   LIB_NETCDF => "/ford1/local/IRIX64/netcdf/lib"   ,
#   INC_NETCDF => "/ford1/local/IRIX64/netcdf/include"  ,
#   LIB_MPI => undef ,
#   INC_MPI => undef ,
#   Arch    => "SGI"   ,             #modfied by arch_
#NCAR blackforest
#   FVData  => "/fs/cgd/csm/inputdata/atm/cam1" ,
#   LMData  => "/fs/cgd/csm/inputdata/lnd/clm2" ,
#   Scratch => "/ptmp" ,
#   LIB_NETCDF => "/usr/local/lib32/r4i4"   ,
#   INC_NETCDF => undef ,
#   LIB_MPI => undef ,
#   INC_MPI => undef ,
#   Arch    => "IBM"   ,             #modfied by arch_
#ORNL Eagle
#   FVData  => "/tmp/gpfs200a/mirin/cam2data/inputdata/atm/cam2" ,
#   LMData  => "/tmp/gpfs200a/mirin/cam2data/inputdata/lnd/clm2" ,
#   Scratch => "/tmp/gpfs200a" ,
#   LIB_NETCDF => "/usr/local/lib32/r4i4"   ,
#   INC_NETCDF => undef ,
#   LIB_MPI => undef ,
#   INC_MPI => undef ,
#   Arch    => "IBM"   ,             #modfied by arch_
#LLNL
# Compaq
#   FVData  => "/nfs/tmp0/mirin/camdata/inputdata/atm/cam1" ,
#   LMData  => "/nfs/tmp0/mirin/camdata/inputdata/lnd/clm2" ,
#   Scratch => "/nfs/tmp0" ,
# IBM Blue/Frost
    FVData  => "/p/gf1/mirin/cam2data/inputdata/atm/cam2" ,
    LMData  => "/p/gf1/mirin/cam2data/inputdata/lnd/clm2" ,
    Scratch => "/p/gf1" ,
    LIB_NETCDF => "/usr/local/netcdf/netcdf-3.5.0/32_64/lib" ,
    INC_NETCDF => undef ,
    LIB_MPI => undef ,
# Compaq
#   INC_MPI => undef ,
#   Arch    => "OSF1"   ,            #modfied by arch_
# IBM Blue/Frost
    INC_MPI => "/usr/lpp/ppe.poe" ,
    Arch    => "IBM"   ,             #modfied by arch_
#NERSC IBM
#   FVData  => "/scratch/scratchdirs/mirin/cam2data/inputdata/atm/cam2" ,
#   LMData  => "/scratch/scratchdirs/mirin/cam2data/inputdata/lnd/clm2" ,
#   Scratch => "/scratch/scratchdirs" ,
#   LIB_NETCDF => "/usr/common/usg/netcdf/3.4/lib" ,
#   INC_NETCDF => "/usr/common/usg/netcdf/3.4/include" ,
#   LIB_MPI => undef ,
#   INC_MPI => "/usr/lpp/ppe.poe" ,
#   Arch    => "IBM"   ,             #modfied by arch_
    OS      => `uname -s | tr -d '\n' ` ,
    Mach    => `uname -m | tr -d '\n' ` ,
    Host    => `uname -n | tr -d '\n' ` ,
    SPMD    => undef ,               #determined by MPI
    SMP     =>  undef,               #determined by OMP
    BaseDir => undef,                #determined by Home
    MODEL_SRCDIR => undef,
    MODEL_EXEDIR => undef,
    MODEL_DATDIR => undef,
    MODEL_DEPDIR => undef,
    MODEL_OBJ    => undef,
    LM_DATDIR    => undef,
#
    PNATS  => 1,
    PCNST  => 5,
    SaveToMSS => 0,
    SaveToSilo => 0,
);
#
# 11, SGI + NAS
# 21, Linux + pgf90
# 31, IBM
my %dep_fields_12 = (Scratch => '/scratch/$LOGNAME',);
#
# Linux + lf95
my %dep_fields_21 = (
   LIB_MPI    => "/usr/local/mpich/lib"     ,
   INC_MPI    => "/usr/local/mpich/include" ,
   LIB_NETCDF => "/usr/local/netcdf/lib"    ,
   INC_NETCDF => "/usr/local/netcdf/include" ,
   Scratch    => '~/scratch',
   Case       => "z18"     ,
   N_SMP     => 2        ,
   PCNST  => 1,
   SaveToMSS => 0,
   SaveToSilo => 0,

);
#
# Linux + pgf90
#
my %dep_fields_22 = (
   LIB_MPI    => "/usr/local/mpich/lib"     ,
   INC_MPI    => "/usr/local/mpich/include" ,
   LIB_NETCDF => "/usr/local/netcdf/lib"    ,
   INC_NETCDF => "/usr/local/netcdf/include" ,
   Scratch    => '~/scratch',
   Case       => "z18"     ,
   N_SMP     => 2        ,  
   SaveToMSS => 0,
   SaveToSilo => 0,
);
#
# Linux + f90
#
my %dep_fields_23 = (
   LIB_MPI     => "/usr/local/mpich-1.2.1-ffc/lib"     ,
   INC_MPI     => "/usr/local/mpich-1.2.1-ffc/include" ,
   LIB_NETCDF  =>  "/usr/local/netcdf-3.4_FFC/lib"     ,
   INC_NETCDF  =>  "/usr/local/netcdf-3.4_FFC/include" ,
   Scratch    => '~/scratch',
   Case       => "z18"     ,
   SaveToMSS => 0,
   SaveToSilo => 0,
);
#
# IBM, phalanx
#
my %dep_fields_31 = (
   LIB_NETCDF  =>  "/usr/local/lib",
   INC_NETCDF  =>  "/usr/local/include",
##   Scratch     => '~/scratch',
##   FVData      => '/g2/home/bwshen//CCM4/data',
   Scratch     => '/scratch',
   FVData      => '/share/fvccm/FVCCM4/data',
   USER_FC     => '',
   N_SMP       => 8        ,
   SaveToMSS => 0,
   SaveToSilo => 0,
);
#
# NERSC IBM
#
my %dep_fields_32 = (
   FVData  => "/u2/mirin/dycore_data/fvccm4" ,
   Scratch => "/scratch" ,
   LIB_NETCDF => "/usr/common/usg/netcdf/3.4/lib" ,
   INC_NETCDF => "/usr/common/usg/netcdf/3.4/include" ,
   LIB_MPI => undef ,
   INC_MPI => "/usr/lpp/ppe.poe" ,
   USER_FC     => '',
   SaveToMSS => 0,
   SaveToSilo => 0,
);

sub new {
    my $that  = shift;
    my $class = ref($that) || $that;
    my $self  = {
        _permitted => \%fields,
        %fields,
    };
    bless $self, $class;
    return $self;
} 
#
#
#
sub print_data_ {

    my $self = shift;
    open (FH, ">conf.log") ;
    print FH "########################################## \n" ;
    print FH "#The parameters you are going to use are:# \n" ;
    print FH "########################################## \n" ;
    foreach $key (sort(keys(%fields))) {
     if ( $self->{$key} eq $fields{$key} ) {
        printf FH "%-15.15s = %s \n",  $key , $self->{$key}  ;
     } else {
      if ( $fields{$key} eq undef ){
      printf FH "%-15.15s = %-40s \n",  $key , $self->{$key};
      }
      else {
      printf FH "%-15.15s = %-40s %15s\n",  $key , $self->{$key}, '<---modified';
      }
     }
    }
    close (FH);
    return ;
}

#
#
#
sub modify {

    my $self = shift;
    my %args;
    if (@_ == 0) {                  # no further arguments
#        print "Null \n";
    }  else {
        %args = @_;                 # use the ones given
    }
    foreach $key (keys(%args)) {
        $self->{$key} = $args{$key} ;
    }
    return ;
}

#
#
#
sub arch_{
my $self = shift;
my ($os) = $self->{'OS'} ;
chomp($os);
my $arch="" ;             #defined here
if    ($os eq IRIX64){
        $arch=SGI ;
}
elsif ($os eq AIX) {
        $arch=IBM ;
}
elsif ($os eq Linux) {
#       my $arch=LINUX ;  #only valid in this {}
        $arch=LINUX ;
}
elsif ($os eq OSF1) {
        $arch=DEC  ;
}
else {
        $arch=CRAY  ;
}

# print "$arch \n";
return $arch;

}

sub change_parm_{
my $self = shift;
my $Arch = $self->{'Arch'};
my $host = $self->{'Host'};
chomp($Arch);
my $USER_FC  = $self->{'USER_FC'};
chomp($USER_FC);

my %tmp;
if ( $Arch eq LINUX ) {
$self->{'LIB_MPI'} = "" ;
$self->{'INC_MPI'} = "" ;
## print "Before: Compiler used $USER_FC \n";
  if ($USER_FC eq pgf90) {
#     print $USER_FC;
     %tmp=(%dep_fields_22);
  } elsif ($USER_FC eq f90) {
#     print "hi $USER_FC";
     %tmp=(%dep_fields_23);
  } elsif ($USER_FC eq lf95) {
#     print $USER_FC;
     %tmp=(%dep_fields_21);
  }
} elsif ($Arch eq SGI) {
#
##  my $host=`uname -n | tr -d '\n' `;
  %tmp=(%dep_fields_12);
} elsif ($Arch eq IBM) {
##  my $hosttmp=`uname -n | tr -d '\n' `;
  $host = substr($host, 0, 4);
  if ($host eq f01n) {
    %tmp=(%dep_fields_31);
    $self->modify(Host => f01n);
  } elsif ($host eq UNDEFINED) {
    %tmp=(%dep_fields_32);
  }
}
else {} ;

$self->modify(%tmp);


#
# read CVS Tag into Ver
#
my $ver = $self->{'Ver'}; 
my $tmp, $tmp1;
if (-e "../CVS/Tag") {
 open (FH, "<../CVS/Tag");
 while (<FH>){
 $buffer=$_;
 if (/ccm/){
    $tmp1=substr($buffer, 6, 5);
    ($tmp) = split (" ", $tmp1);   #remove trailing spaces, 11_62 or 12_6
    $ver="fvcam_${tmp}";
    print "$ver\n" ;
 }
 }
$self->modify(Ver => $ver);
}



}
#
#
#
sub read_conf{

use FileHandle;
my $FH = new FileHandle;

my $self = shift;
my ($conf_log) = @_;
if ($conf_log eq "-" ){
  $FH=STDIN;
} else{
  open ($FH, "<$conf_log");
}
#print "FH $FH \n";
while ($buffer = <$FH>){
 next if $buffer =~ /^#/ ;
 @pairs = split(/</, $buffer); 
 ($keytmp, $valuetmp) = split (/=/, $pairs[0]); 
 ($key) = split (" ", $keytmp);      #remove the trailing spaces
 ($value) = split (" ", $valuetmp);  #remove the leading spaces
 chomp($key); chomp($value);
# print "$key . $value \n";
 if ( $value ne "" ) {
    $self->{$key}=$value;
 }
# print "$key .  $self->{$key} \n" ;
}
 return;
}

#
#
#
sub user_input{
  my $self = shift;

  if ($self->{'Arch'} == LINUX) {
  print "Compiler ($self->{'USER_FC'}):";
  chomp($ans = <STDIN> );
  if ($ans ne "" ) {  ($self->{'USER_FC'}) = split(" ", $ans) };
  $self->change_parm_;
  }

  print "Simulation Dynamics: Options are  1:LR   3:EUL ($self->{'Dyn'}):" ;
  chomp($ans = <STDIN> ); $ans and ($self->{'Dyn'}) = split(" ", $ans);

  print "Simulation LAND Model: Options are  1:clm   2:lsm ($self->{'Phy'}):" ;
  chomp($ans = <STDIN> ); $ans and ($self->{'Phy'}) = split(" ", $ans);

  print "Ver  id ($self->{'Ver'}):" ;
  chomp($ans = <STDIN> ); $ans and ($self->{'Ver'}) = split(" ", $ans);

  print "Case id ($self->{'Case'}):";
  chomp($ans = <STDIN> ); $ans and ($self->{'Case'}) = split(" ", $ans);
#  print "case $self->{'Case'} \n";


  print "FVData ($self->{'FVData'}):";
  chomp($ans = <STDIN> ); $ans and ($self->{'FVData'}) = split(" ", $ans);


  print "SCRATCH ($self->{'Scratch'}):";
  chomp($ans = <STDIN> ); $ans and ($self->{'Scratch'}) = split(" ", $ans);

  if ($self->{'Arch'} ne LINUX) {
  print "LIB_NETCDF ($self->{'LIB_NETCDF'}):";
  chomp($ans = <STDIN> ); $ans and ($self->{'LIB_NETCDF'}) = split(" ", $ans);

  print "INC_NETCDF ($self->{'INC_NETCDF'}):";
  chomp($ans = <STDIN> ); $ans and ($self->{'INC_NETCDF'}) = split(" ", $ans);
  }

  print "OpenMP on? ($self->{'OMP'}):" ;
#  chomp($ans = <STDIN> ); $ans and ($self->{'OMP'}) = split(" ", $ans);
#                          ^^^^ $ans must be an non-zero number!
  chomp($ans = <STDIN> ); 
  if ($ans ne "" ) {  ($self->{'OMP'}) = split(" ", $ans) };
  if ($self->{'OMP'} eq 1 ) {
     print "Numbers of OpenMP threads ($self->{'N_SMP'}):";
     chomp($ans = <STDIN> ); $ans and ($self->{'N_SMP'}) = split(" ", $ans) ;
  } else {
    $self->{'N_SMP'} = 1;
  }
  

  $self->{'SMP'}=$bool{$self->{'OMP'}};

  print "MPI on? ($self->{'MPI'}): " ;
  chomp($ans = <STDIN> ); 
  if ($ans ne "" ) { ($self->{'MPI'}) = split(" ", $ans) };
  if ($self->{'MPI'} eq 1 ) {
     print "Numbers of MPI processes ($self->{'N_MPI'}):";
     chomp($ans = <STDIN> ); $ans and ($self->{'N_MPI'}) = split(" ", $ans);
     print "Use oned or twod decomposition ($self->{'Decomp'}):";
     chomp($ans = <STDIN> ); $ans and ($self->{'Decomp'}) = split(" ", $ans);
     if ($self->{'Decomp'} eq "twod" ) {
       print "Numbers of MPI processes ($self->{'NPR_YZ'}):";
       chomp($ans = <STDIN> ); $ans and ($self->{'NPR_YZ'}) = split(" ", $ans);
     }
       
     if ($self->{'Arch'} ne LINUX) {
     print "LIB_MPI ($self->{'LIB_MPI'}):";
     chomp($ans = <STDIN> ); $ans and ($self->{'LIB_MPI'}) = split(" ", $ans);

     print "INC_MPI ($self->{'INC_MPI'}):";
     chomp($ans = <STDIN> ); $ans and ($self->{'INC_MPI'}) = split(" ", $ans);
     }
  } else {
    $self->{'N_MPI'} = 1;
  }

  $self->{'SPMD'}=$bool{$self->{'MPI'}};


  print "PNATS? ($self->{'PNATS'}): " ;
  chomp($ans = <STDIN> );
  if ($ans ne "" ) { ($self->{'PNATS'}) = split(" ", $ans) };


  print "PCNST? ($self->{'PCNST'}): " ;
  chomp($ans = <STDIN> );
  if ($ans ne "" ) { 
     if ($ans > 1 ){$ans=5};
    ($self->{'PCNST'}) = split(" ", $ans) 
  };


 



  

}

#
#
#
sub model_vars_{
my $self = shift;
my $ver  = $self->{'Ver'};
my $case = $self->{'Case'};
my $obj  = $self->{'Obj'};
my $basedir = dirname($self->{'Home'});
$self->{'MODEL_DATDIR'} = $self->{'FVData'} ;
$self->{'LM_DATDIR'} = $self->{'LMData'} ;
$self->{'MODEL_SRCDIR'} = "${basedir}/src" ;
$self->{'MODEL_EXEDIR'} = "${basedir}/${ver}/${case}" ;
$self->{'MODEL_DEPDIR'} = "${basedir}/tools/makdep" ;
$self->{'MODEL_OBJ'}    = "${basedir}/${ver}/${case}/${obj}" ;

}
#
#
#
sub mkpath_ {       #can not be mkpath
use File::Path;
my $self = shift;

mkpath ([$self->{'MODEL_EXEDIR'}, $self->{'MODEL_OBJ'}], 1, 0711) ;


}

#
#
#
sub headerFiles_{
my $self = shift;

chdir $self->{'MODEL_OBJ'};

my $dyn  = $self->{'Dyn'}  ;
my $arch = $self->{'Arch'} ;
my $mpi  = $self->{'MPI'}  ;
my $case = $self->{'Case'} ;
my $pnats = $self->{'PNATS'};
my $pcnst = $self->{'PCNST'};
my $decomp = $self->{'Decomp'};

local $idim, $jdim, $kdim;
if ($dyn eq 1 ) {
   setdim($case) ; 
} elsif ($dyn eq 3) {
#   setdim_eul($case) ;
    $idim=128;
    $jdim=64;
} else {}

misc_h_($dyn, $arch, $mpi, $decomp);
params_h_($dyn, $idim, $jdim, $kdim, $pnats, $pcnst, $decomp);
preproc_h_($dyn, $idim, $jdim, $kdim);

chdir $self->{'Home'};


}

#
#
#
sub setdim() {
my ($case) = @_ ;

#@levels = (96, 55, 32, 30, 18 ); # supported levels
 @levels = (30, 26, 18); # supported levels

my $hres = substr($case,0,1);
my $vres = substr($case,1,2);

#my $idim, jdim;

# Figure out horizontal resolution
# --------------------------------
if ( lc "$hres" eq "a" ) {
     $idim = 72;
     $jdim = 46;
   } elsif ( lc "$hres" eq "b" ) {
     $idim = 144;
     $jdim = 91;
   } elsif ( lc "$hres" eq "c" ) {
     $idim = 288;
     $jdim = 181;
   } elsif ( lc "$hres" eq "d" ) {
     $idim = 576;
     $jdim = 361;
   } elsif ( lc "$hres" eq "p" ) {
     $idim = 100;
     $jdim = 61;
   } elsif ( lc "$hres" eq "z" ) {
     $idim = 24;
     $jdim = 19;
   } else {
     die "Invalid FVCAM resolution $res";
   }
#



# Figure out vertical resolution
# ------------------------------
   my $ok = 0;
   my $lev;
   foreach $lev ( @levels ) {
     if ( $vres == $lev ) {
          $ok = 1;
          $kdim = $vres;
        }
   }
   die "Unsupported number of FVCAM levels $vres" unless ( $ok );


}


#
#
#
sub misc_h_() {
my ($dyn, $arch, $mpi, $decomp) = @_ ;
if ($arch eq DEC) {$arch=OSF1};
if ($arch eq IBM) {$arch=AIX};
open (FH, ">misc.h");
#
if ($dyn eq  11111 ) {
print FH <<"EOF"; 
#ifndef MISC_SET
#define MISC_SET
#define $arch
#undef  PVP 
#define REALTYPE MPI_DOUBLE_PRECISION
#undef  COUP_SOM
#undef  COUP_CSM
#define NCPREC   NF_FLOAT
#define FORTFFT
#define STAGGERED
EOF
}
#
elsif ($dyn eq 1) {
print FH <<"EOF"; 
#ifndef MISC_SET
#define MISC_SET
#undef  PVP
#undef  COUP_SOM
#undef  COUP_CSM
#define FORTFFT
#define STAGGERED
#undef  PERGRO
#define PCWDETRAIN
EOF
}
#
elsif ($dyn eq 33333 ) {
print FH <<"EOF"; 
#ifndef MISC_SET
#define MISC_SET
#define $arch
#undef  COUP_CSM
#define FORTFFT
EOF
}
#
elsif ($dyn eq 3 ) {
print FH <<"EOF"; 
#ifndef MISC_SET
#define MISC_SET
#undef  COUP_CSM
#define FORTFFT
#undef  PERGRO
#define PCWDETRAIN
EOF
}
#
else { } ;

if ($mpi == 0 ){
 print FH "#undef  SPMD \n" ;
}
else {
 print FH "#define SPMD \n" ;
}

if ($arch eq SGI) {
  print FH "#define SGI_FFT \n";
}

print FH "#endif \n" ;


close(FH);
return ;
}

#
#
#

sub params_h_ {
my ($dyn, $idim, $jdim, $kdim, $pnats, $pcnst, $decomp) = @_ ;
open (FH, ">params.h");
#
if ( $dyn eq 1 ) {
print FH <<"EOF";
#ifndef PARAMS_SET
#define PARAMS_SET
#define PNATS  $pnats
#define PCNST  $pcnst
#define PLEV   $kdim
#define PLEVR  $kdim
#define PLON   $idim
#define PLAT   $jdim
#define PTRM  42 
#define PTRN  42
#define PTRK  42
#define PCOLS 16
#endif
EOF
}
elsif ( $dyn eq 3 ) {
print FH <<"EOF";
#ifndef PARAMS_SET
#define PARAMS_SET
#define PCNST  1
#define PNATS  1
#define PLEV   30
#define PLEVR  30
#define PLON   $idim
#define PLAT   $jdim
#define PTRM   42
#define PTRN   42
#define PTRK   42
#endif
EOF
}
else {};

if ($decomp eq "twod" ){
 print FH "#define TWOD_YZ \n" ;
}

close(FH);

return;

}

#
#
#
sub preproc_h_ {
my ($dyn, $idim, $jdim, $kdim) = @_ ;
open (FH, ">preproc.h");
#
if ( $dyn eq 1 ) {
print FH  <<"EOF" ;
#ifndef PREPROC_SET
#define PREPROC_SET
#define COUP_CAM
#define LSMLON $idim
#define LSMLAT $jdim
#endif
EOF
}
#
elsif ( $dyn eq 3 ) {
print FH  <<"EOF" ;
#ifndef PREPROC_SET
#define PREPROC_SET
#define COUP_CAM
#define LSMLON $idim
#define LSMLAT $jdim
#endif
EOF
}
#
else {} ;
close(FH);
return ;
}

#
#
#
sub makefile_n_filepath_{

my $self = shift;

chdir $self->{'MODEL_OBJ'};
my $dyn = $self->{'Dyn'};
my $phy = $self->{'Phy'};
my $model_srcdir = $self->{'MODEL_SRCDIR'};
my $omp = $self->{'OMP'};
my $source = $self->{'BaseDir'};
my $home = $self->{'Home'};
my $os  =$self->{'OS'};
if ($os eq AIX) {
  $os=rs6000_sp;
}

Filepath_($dyn, $phy, $model_srcdir, $os);
$self->Make_macros_();
Makefile_($omp, $source, $home);

chdir $self->{'MODEL_EXEDIR'};
$self->mMakefile_();

chdir $self->{'Home'};
}



#
sub Filepath_ {
my ($dyn, $phy, $MODEL_SRCDIR, $os) = @_ ;
open (FH, ">Filepath");
#
if ( $dyn eq 1 ){
   if ( $phy eq 1 ){
print FH  <<"EOF" ;
$MODEL_SRCDIR/mods
$MODEL_SRCDIR/dynamics/fv
$MODEL_SRCDIR/control
$MODEL_SRCDIR/physics/cam1
$MODEL_SRCDIR/../../../atmlnd_share
$MODEL_SRCDIR/../../../csm_share
$MODEL_SRCDIR/ocnsice/dom
$MODEL_SRCDIR/utils
$MODEL_SRCDIR/../../../utils/pilgrim
$MODEL_SRCDIR/../../../utils/timing
$MODEL_SRCDIR/../../../lnd/clm2/src/main
$MODEL_SRCDIR/../../../lnd/clm2/src/biogeophys
$MODEL_SRCDIR/../../../lnd/clm2/src/biogeochem
$MODEL_SRCDIR/../../../lnd/clm2/src/mksrfdata
$MODEL_SRCDIR/../../../lnd/clm2/src/ecosysdyn
$MODEL_SRCDIR/../../../lnd/clm2/src/riverroute
$MODEL_SRCDIR/../../../ice/csim4
EOF
  }
  elsif ($phy eq 2) {
print FH  <<"EOF" ;
$MODEL_SRCDIR/mods
$MODEL_SRCDIR/../../utils/pilgrim
$MODEL_SRCDIR/../../utils/timing
$MODEL_SRCDIR/dynamics/fv
$MODEL_SRCDIR/control
$MODEL_SRCDIR/ocnsice/dom
$MODEL_SRCDIR/physics/cam
$MODEL_SRCDIR/../../csm_share
$MODEL_SRCDIR/../../atmlnd_share
$MODEL_SRCDIR/../../lnd/lsm/src
$MODEL_SRCDIR/utils
EOF
  }
  else {};
}
#
elsif ($dyn eq 3 ){
print FH <<"EOF" ;
$MODEL_SRCDIR/../../utils/timing
$MODEL_SRCDIR/mods4eul
$MODEL_SRCDIR/control
$MODEL_SRCDIR/dynamics
$MODEL_SRCDIR/dynamics/eul
$MODEL_SRCDIR/ocnsice/dom
$MODEL_SRCDIR/physics/cam
$MODEL_SRCDIR/../../csm_share
$MODEL_SRCDIR/atmlnd_share
$MODEL_SRCDIR/../../ldn/lsm
$MODEL_SRCDIR/utils
EOF
}
#
else {} ;
close(FH);
return ;
}

#
#
#


sub Make_macros_ {
my $self = shift;
open (FH, ">Make.macros");
print FH   <<"EOF" ;
SPMD:=$self->{'SPMD'}
SMP:=$self->{'SMP'}
LIB_NETCDF:=$self->{'LIB_NETCDF'}
INC_NETCDF:=$self->{'INC_NETCDF'}
LIB_MPI:=$self->{'LIB_MPI'}
INC_MPI:=$self->{'INC_MPI'}
MODEL_EXEDIR:=$self->{'MODEL_EXEDIR'}
EXENAME:=$self->{'EXE'}
USER_FC:=$self->{'USER_FC'}
#LIB_ESMF:=$self->{'MODEL_SRCDIR'}"/utils/esmf/lib/libO"
EOF
#
close(FH);
#
open (FH, ">Rootdir");
print FH "$self->{'MODEL_SRCDIR'}/../../../../\n";
close(FH);

return;

}

#
#
#
sub Makefile_ {
my ($omp, $source, $Home) = @_ ;
if ($omp eq 1 ) {
  $mp="-mp";
  $openmp="--openmp";
} else {
 $mp="" ;
 $openmp="";
};
my $makefile_output="Makefile";
my $makefile="${Home}/Makefile" ;
print "Makefile is $makefile \n" ;
print "Source is $source \n" ;
open (FH, ">$makefile_output");
print FH "include Make.macros \n";
#1close(FH);
#1
#1`cat ${Home}/Makefile >>  Makefile`

if ( ! -e "${makefile}") {
if (-e "${source}/bld/Makefile"){
#   $makefile = "${source}/bld/Makefile.tmpl";  #bug?
   $makefile = "${source}/bld/Makefile";
 }
elsif (-e "${Home}/Makefile"){
   $makefile = "${Home}/Makefile";
 }
else {die "ERROR: ${Home}/Makefile not found \n";}
}
print "Makefile is $makefile \n" ;
open (IN, "<${makefile}");
my $fc;
while (<IN>) {
   @fld=split(' ', $_, 9999);
   $fc=substr(@fld[0], 0, 1);
   if ($fc eq "#") {
      print FH $_  ;
   } else{
      s/-mp/$mp/o;
      s/--openmp/$openmp/o;     #Lahey Fortran
      s/pgcc/gcc -DUSE_GCC -DDISABLE_TIMERS/o;   #for DAO linux machines only
      if ($omp eq 0 ) {
#        s/cc/cc -DUSE_GCC -DDISABLE_TIMERS/o;#disable OpenMP threads, not good
      }
      if (/^F90FLAGS/) {
        chomp();
        print FH "$_ -I\$\(INC_NETCDF\) \n";
###      } elsif (/^BASE_FFLAGS/){
###        chomp();
###        print FH "$_ -Mextend \n";
      } elsif (/^CFLAGS/){
        s/-fast//o;                       #LINUX only
###     s/-DLINUX -fast/-DLINUX/o;        #OK! better?
        (print FH  $_);
#        chomp() ;                                 #11/27/2001
#        (print FH  "$_ -DCOUP_CAM \n");           #11/27/2001
      }
        elsif (/^FFLAGS/){
        s/-DIRIX64/-DIRIX64 -DSGI -OPT:Olimit=4000/o;
        (print FH  $_);
#        chomp() ;                                 #11/27/2001
#        (print FH  "$_ -DCOUP_CAM \n ");           #11/27/2001
      }  elsif (/^LDFLAGS/){
        s/-64/-64 -lcomplib.sgimath -lfastm -s/o;   #SGI only
        (print FH  $_);
      }
      else{
      (print FH  $_);
      }
   };
}
close(FH);
close(IN);
}
#
#
#
sub namelist_{

#chdir $MODEL_EXEDIR;
}

#
#
#
sub run_scripts_{

my $self = shift;

chdir $self->{'MODEL_EXEDIR'};

my $Home = $self->{'Home'};

-e ("$Home/script.m4") or die "ERROR: $Home/script.m4 not found \n" ;
my $arch    = $self->{'Arch'};
my $mach    = $self->{'Mach'};
my $decomp  = $self->{'Decomp'};
my $mpi     = $self->{'MPI'};
my $n_mpi   = $self->{'N_MPI'};
my $npr_yz  = $self->{'NPR_YZ'};
my $omp     = $self->{'OMP'};
my $n_smp   = $self->{'N_SMP'};
my $basedir = $self->{'BaseDir'};
my $ver     = $self->{'Ver'};
my $n_total = $self->{'N_TOTAL'};
my $case    = $self->{'Case'};
my $exe     = $self->{'EXE'};
my $fvdata  = $self->{'FVData'};
my $scratch = $self->{'Scratch'};
my $dyn = $self->{'Dyn'};
my $host = $self->{'Host'};
my $pcnst = $self->{'PCNST'};
my $savetomss = $self->{'SaveToMSS'};
my $savetosilo = $self->{'SaveToSilo'};
my $lmdata  = $self->{'LMData'} ;



if ("$mpi" eq 1 ){
   if ("$omp" eq 1 ) {
       $n_total = $n_mpi * $n_smp ;
   } else {
       $n_total = $n_mpi;
   }
} elsif ("$omp" eq 1 ){
        $n_total = $n_smp ;
} else  {$n_total = 1; }
$self->{'N_TOTAL'} = $n_total;




my $m4     = "define.m4";
my $m4_tmpl= "$m4.tmpl" ;
define_m4_tmpl($m4_tmpl, $mpi, $n_mpi, $npr_yz, $omp, $n_smp, $basedir, $ver, $case, $exe, $fvdata, $scratch, $dyn, $host, $pcnst, $savetomss, $savetosilo, $lmdata, $decomp, $n_total) ;
@SCRIPTS=("GENERIC", "WEEKLY", "MONTHLY");
foreach my $type (@SCRIPTS){
  if ( "$type" eq  "GENERIC" ) {
                  $script="fvcam.j";
  } elsif ( "$type" eq "WEEKLY" ) {
                  $script="weekly.j.tmpl" ;
  } elsif ( "$type" eq "MONTHLY" ) {
                  $script="monthly.j.tmpl" ;
  } else  {
                  $script="tmp" ;;
  }

  sed_m4($arch, $mach, $type, $m4_tmpl, $m4) ;

  `/usr/bin/m4 -B1000000 $Home/script.m4 > $script` ;

  chmod 0755, $script ;

}

unlink ($m4);
unlink ($m4_tmpl);

chdir $self->{'Home'};
return;
}
#
#
#

sub define_m4_tmpl {
 my ($tmpl, $mpi, $n_mpi, $npr_yz, $omp, $n_smp, $basedir, $ver, $case, $exe, $fvdata, $scratch, $dyn, $host, $pcnst, $savetomss, $savetosilo, $lmdata, $decomp, $n_total) = @_ ;

my $hres = substr($case,0,1);
my $vres = substr($case,1,2);
my $fv_ncdata, $fv_bndtvs ;
my $fv_dtime, $fv_lsmgeo;
my $fv_clm_surf;

my $trace_gas='.true.' ;
if ($pcnst eq 1 ){
   $trace_gas='.false.';
}

my $n_node=$n_mpi;   #IBM only 
if ($n_node eq 0 ){
  $n_node = 1;
} 

# --------------------------------
if ($dyn eq 1 ) {
#$fv_dtime = 3600;
$fv_dtime = 1800;
$fv_lsmgeo = 33;
if ( lc "$hres" eq "a" ) {
     $fv_ncdata="fv.4x5.nc";
     $fv_bndtvs="sst_HadOIBl_bc_4x5_clim_c020411.nc";
     $fv_clm_surf="clms_4x5_c020412.nc";
     if ( lc "$vres" eq "26") {
        $fv_ncdata="cami_0000-09-01_4x5L26_c020430.nc" ;
     }
   } elsif ( lc "$hres" eq "b" ) {
     $fv_ncdata="fv.2x2.5.nc";
     $fv_bndtvs="sst_HadOIBl_bc_2x2.5_clim_c020531.nc";
     $fv_clm_surf="clms_2x2.5_c021004_USGSsm.nc";
     if ( lc "$vres" eq "26") {
        $fv_ncdata="cami_0000-09-01_2x2.5_L26_c021209_USGS.nc";
     }
   } elsif ( lc "$hres" eq "c" ) {
     $fv_ncdata="fv.1x1.25.nc";
     $fv_bndtvs="sst_HadOIBl_bc_1x1.25_clim_c021210.nc";
     $fv_clm_surf="clms_1x1.25_Navy_c021210.nc";
     if ( lc "$vres" eq "26") {
        $fv_ncdata="cami_0000-09-01_1x1.25_L26_c021127_Navy.nc";
     }
   } elsif ( lc "$hres" eq "d" ) {
     $fv_ncdata="fv.0.5x0.625.nc";
     $fv_bndtvs="amipbc_sst_sic_361x576_clim.nc";
     $fv_clm_surf="surface-data.576x361.nc";
     if ( lc "$vres" eq "26") {
        $fv_ncdata="cami_0000-09-01_0.5x0.625_oro_L26_c020430.nc";
     }
   } elsif ( lc "$hres" eq "p" ) {
     $fv_ncdata="fv.3x3.6.nc";
     $fv_bndtvs="lr.SST3x3.6.nc";
     $fv_clm_surf="not_found";
   } elsif ( lc "$hres" eq "z" ) {
     $fv_ncdata="fv.10x15.nc";
     $fv_bndtvs="sst_HadOIBl_bc10x15_clim.c020411.nc";
     $fv_clm_surf="clms_10x15_c020412.nc";
     if ( lc "$vres" eq "26") {
        $fv_ncdata="cami_0000-09-01_10x15_L26_c020430.nc" ;
     }
   } else {
     die "Invalid FVCAM resolution $res";
   }
} elsif ($dyn eq 3 ) {
 $fv_dtime = 1200;
 $fv_lsmgeo = 23;
 $fv_ncdata="SEP1.T42L30.051700.definesurf.nc";
 $fv_bndtvs="T42M5079.nc";
}
else {}


#    print "tmpl ", $tmpl, "\n";
    open (FH, ">$tmpl") or die "Cannot open $tmpl" ;
#left alignment
    print FH <<"EOF" ;
dnl
define([%FV_ARCH%], [%arch%])dnl
dnl
define([%FV_MACH%], [%mach%])dnl
dnl
define([%FV_HOST%], [%$host%])dnl
dnl
define([%FV_TYPE%], [%type%])dnl
dnl
define([%FV_STAGE%], [%staging%])dnl
dnl
define([%FV_DECOMP%], [%$decomp%])dnl
dnl
define([%FV_MPI%], [%$mpi%])dnl
dnl
define([%FV_N_MPI%], [%$n_mpi%])dnl
dnl
define([%FV_NPR_YZ%], [%$npr_yz%])dnl
dnl
define([%FV_NODE%], [%$n_node%])dnl
dnl
define([%FV_OMP%], [%$omp%])dnl
dnl
define([%FV_N_SMP%], [%$n_smp%])dnl
dnl
define([%FV_BASEDIR%], [%$basedir%])dnl
dnl
define([%FV_VER%], [%$ver%])dnl
dnl
define([%FV_CASE%], [%$case%])dnl
dnl
define([%FV_EXE%], [%$exe%])dnl
dnl
define([%FV_DATA%], [%$fvdata%])dnl
dnl
define([%FV_NCDATA%], [%$fv_ncdata%])dnl
dnl
define([%FV_BNDTVS%], [%$fv_bndtvs%])dnl
dnl
define([%FV_SCRATCH%], [%$scratch%])dnl
dnl
define([%FV_DTIME%], [%$fv_dtime%])dnl
dnl
define([%FV_LSMGEO%], [%$fv_lsmgeo%])dnl
dnl
define([%FV_TRACE_GAS%], [%$trace_gas%])dnl
dnl
define([%FV_SaveToMSS%], [%$savetomss%])dnl
dnl
define([%FV_SaveToSilo%], [%$savetosilo%])dnl
dnl
define([%FV_LMDATA%], [%$lmdata%])dnl
dnl
define([%FV_CLM_SURF%], [%$fv_clm_surf%])dnl
dnl
define([%NCPUGEN%], [%$n_total%])dnl
dnl --------------------------------------
dnl --------------------------------------
EOF
    close (FH);

    return;
}
#
#
#
sub sed_m4 {
    my ($arch, $mach, $type, $old, $new) = @_ ;

#    print "type  ", "$type ", " \n" ;

    open(OLD, "< $old");
    open(NEW, "> $new");

    while (<OLD>) {
        s/arch/$arch/o;
        s/mach/$mach/o;
        s/type/$type/o;
        (print NEW  $_);
    }

    close(OLD);
    close(NEW);

}

#
#
#
sub mMakefile_{

my $self = shift;
my ($arch) =($self->{'Arch'});
my ($os)  = ($self->{'OS'})  ;
my ($exedir) = ($self->{'MODEL_EXEDIR'});
my ($model_srcdir)  = ($self->{'MODEL_SRCDIR'});
my $comm="gmake -j1";

if    ($arch eq SGI){
      $comm='gmake -j8; gmake -j8';
}
elsif ($arch eq IBM) {
      $comm='gmake -j4; gmake -j4';
}
elsif ($arch eq LINUX) {
      $comm='gmake -j2';
}
elsif ($arch eq DEC) {
#      $comm='gmake -j2; gmake -j2';
}
else {
#      $comm='gmake -j2; gmake -j2';
}

if ($os eq AIX) {             #phalanx
    $os=rs6000_sp;
}

open (FH1, ">Makefile");
#print FH1 "MF_ARCH :=$os \n";
#print FH1 "MF_DIR  :=$model_srcdir/utils/esmf\n";
#print FH1 "LIB_ESMF:=$model_srcdir/utils/esmf/lib/libO\n" ;
print FH1 "ALL: \n" ;
#print FH1 "\t (export MF_ARCH=\$(MF_ARCH) MF_DIR=\$(MF_DIR);  cd \$(MF_DIR); gmake BOPT=O)\n";
print FH1 "\t (cd obj; $comm )\n";
print FH1 "clean: \n";
print FH1 "\t (cd obj; gmake clean) \n";
print FH1 "cleanall: \n";
print FH1 "\t (cd obj; gmake clean; rm -rf *.d) \n";
close(FH1);

open (FH1, ">gmake.j");
print FH1 <<"EOF";
#!/bin/csh -f
# ------------------------------
#PBS -l ncpus=8
#PBS -l walltime=1:00:00
#PBS -l mem=4gb
#PBS -S /bin/csh
#PBS -V
#PBS -j eo
# ------------------------------
cd $exedir
gmake
EOF
close(FH1);
chmod 0755, "gmake.j" ;

}


#
#
#


sub usage_{


my $self = shift;
my ($ver, $case) = ($self->{'Ver'}, $self->{'Case'});
my ($arch) =($self->{'Arch'});
my $host = $self->{'Host'}; 
my $run_opt="qsub fvcam.j" ;
my $gflags="";

if    ($arch eq SGI){
       $run_opt="fvcam.j";
       $gflags="-j4 -e";
}
elsif ($arch eq IBM) {
       $run_opt="llsubmit fvcam.j";
       $gflags="-j2 -e";
}
elsif ($arch eq LINUX) {
       $run_opt="fvcam.j";
}
elsif ($arch eq DEC) {
#        $arch=DEC  ;
}
else {
#        $arch=CRAY  ;
}



print  << "EOF";

======Quick Start ========================================================
cd ../$ver/$case
gmake 
$run_opt 
==========================================================================

EOF

}


#
#
#

sub makdep_tool_{
#-----------------------------------------------------------------------
#
# cd to $MODDEL_DEPDIR (default setting is "../tools/makdep) to create
# an alternative dependency generator to cpp.  Running gmake will create
# an executable named "makdep" in $(HOME)/bin.  This directory will need
# to be in one's $PATH in order to build cam with the standard Makefile.
#
#-----------------------------------------------------------------------
my $self = shift;
my ($model_depdir, $Home) = ($self->{'MODEL_DEPDIR'}, $self->{'Home'});
use File::Copy ;

-e ("$model_depdir") or die "ERROR:  $model_depdir not found\n";

my $tmp="${Home}/tmp";
mkdir ($tmp, 0711) ;
#copy ("${model_depdir}/*" , $tmp);
`cp -r ${model_depdir}/*  $tmp` ;
chdir($tmp) ;

chdir($Home)         ;

`rm -rf $tmp `       ;


}




#
#
#


package main;

use  Cwd;
use File::Basename;
use File::Path ;
use Getopt::Std;
use Env;        #for $USER

sub mss_shell {
my ($src, $exedir) = @_;

my $tmpdir="$src/mods";
print "tmpdir $tmpdir \n";
mkdir ($tmpdir, 0711) ;
chdir $tmpdir;


    open(OLD, "< ../../src/control/ioFileMod.F90");
    open(NEW, "> ioFileMod.F90");
    while (<OLD>) {
        s/mswrite/${exedir}\/mswrite/o;
        s/msread/${exedir}\/msread/o;
        s/system_cmd/system_cmd_f/o;
        (print NEW  $_);
    }
    close(OLD);
    close(NEW);


    open(OLD, "< ../../../../../models/lnd/clm2/src/main/fileutils.F90");
    open(NEW, "> fileutils.F90");
    while (<OLD>) {
        s/mswrite/${exedir}\/mswrite/o;
        s/msread/${exedir}\/msread/o;
        s/system_cmd/system_cmd_f/o;
#        s/100a/512a/o;
        (print NEW  $_);
    }
    close(OLD);
    close(NEW);

    open(NEW, "> system_cmd_f.F90");
    print NEW " function system_cmd_f(text) \n";
    print NEW " character(len=\*), intent(in) :: text \n";
    print NEW " character\*256                :: buf  \n";
    print NEW " buf=trim(text)   \n";
    print NEW " system_cmd_f=system_cmd(buf) \n" ;
    print NEW " return \n";
    print NEW " end function system_cmd_f \n";

    close (NEW);



}


getopts('if:');

my $interactive;
my $conf_log;
#if ($#ARGV ge -1) 
#{
#  $interactive=0 unless ($interactive =$opt_i);
#}

$interactive=0 unless ($interactive =$opt_i);
$conf_log   ="conf.log" unless ($conf_log =$opt_f);
## print "interactive $interactive \n";
 print "conf_log $conf_log \n";

# $ccm = CCM_Env->new->Dyn(5);           
$ccm = CCM_Env->new;            #wrong
#$ccm = CCM_Env->new;                     #works, instance method
#$ccm = new CCM_Env;                     #works, class method

#   print "Dyn is    $ccm->Dyn    \n" ;      #wrong
# print "Dyn is ", $ccm->{'Dyn'}, " \n" ;      #correct

#
#
#
$ccm->{'BaseDir'} = dirname ($ccm->{'Home'}) ;

# check Arch
#
$arch=$ccm->arch_ ;
chomp($arch);
$ccm->modify(Arch=> $arch);

#
#
#
### if ( $arch eq LINUX ) {
###   print "OMP is set to be 0. \n" if ($ccm->{'OMP'} eq 1);
### #  $ccm->{'OMP'} = 0;
###   $ccm->modify(OMP=> 0);
###   $ccm->modify(Scratch => '~/scratch');
### #  $ccm->modify(FVData  => '/home/bwshen/CCM4/data/FVCCM4');
###   $ccm->modify(N_SMP   => 1);
### }

# print  "$ccm->{'OMP'} \n" ;


if ( $arch eq LINUX ) {
  $ccm->{'USER_FC'}=lf95;
} elsif ( $arch eq SGI ) {
  $ccm->{'USER_FC'}=f90;
}
else {}

#
# set LIB and INC
#
$ccm->change_parm_;   #could be modified in  user_input()

#
# determine the value of SPMD according to the value of MPI
#
$ccm->{'SPMD'}=$bool{$ccm->{'MPI'}}; #could be modified in user_input()

# interactive
if ($interactive eq 1 ) {
  if (-e $conf_log) {
   $ccm->read_conf($conf_log);
  } else {
  print "Warning:################### \n";
  print "Configuration file $conf_log not found \n";
  print "########################### \n";
  }
#  print "MPI $ccm->{'MPI'} \n" ;
  $ccm->user_input;
} elsif ( $interactive eq 0 ){
  if (-e  $conf_log || $conf_log eq '-') {
        $ccm->read_conf($conf_log);
        $ccm->{'SPMD'}=$bool{$ccm->{'MPI'}}; # 07/10/2001
  } 
} else {}



#
# define MODEL variables
#
$ccm->model_vars_;

#
# For a MPI mode,
#
if ($ccm->{'MPI'} eq 1 ){
  mss_shell($ccm->{'MODEL_SRCDIR'}, $ccm->{'MODEL_EXEDIR'});
}
#
#
#
$ccm->mkpath_;

#
$ccm->headerFiles_;
$ccm->makefile_n_filepath_;
#$ccm->namelist_;
$ccm->run_scripts_;


#
#
#
$ccm->print_data_;
#`/bin/more conf.log` ;
	

$ccm->usage_;

$ccm->makdep_tool_;



