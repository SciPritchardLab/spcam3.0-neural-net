divert(-1)dnl

dnl ---------------------------------------------------
dnl  m4 source for generating generic job script and
dnl  other special purpose scripts such as for monthly
dnl  and weekly jobs. 
dnl ---------------------------------------------------

changequote([%,%])dnl

dnl ---------------------------------------
dnl  Suffix for generic, monthly or weekly
dnl ---------------------------------------

define([%FV_SUFFIX_SAVE%], [%dnl
ifelse(FV_TYPE, MONTHLY, [%.${mmsave}%], [%dnl
ifelse(FV_TYPE, WEEKLY,  [%.${mmsave}-${wksave}%], [%dnl
ifelse(FV_TYPE, GENERIC, [%%], [%dnl
ifelse(FV_TYPE, DYCORE, [%%])dnl
%])dnl
%])dnl
%])dnl
%])dnl

define([%FV_DAYS%], [%dnl
ifelse(FV_TYPE, MONTHLY, [%${days}%], [%dnl
ifelse(FV_TYPE, WEEKLY,  [%${days}%], [%dnl
ifelse(FV_TYPE, GENERIC,  [%1%], [%dnl
ifelse(FV_TYPE, DYCORE, [%1%])dnl
%])dnl
%])dnl
%])dnl
%])dnl

dnl --------------------
dnl Define filenames for IC
dnl --------------------
define([%FV_SSTDATA%], [%dnl
ifelse(FV_TYPE, MONTHLY, [%'${s_name}_${s_res}.${s_freq}.${s_date}.nc'%], [%dnl
ifelse(FV_TYPE, WEEKLY,  [%'${s_name}_${s_res}.${s_freq}.${s_date}.nc'%], [%dnl
ifelse(FV_TYPE, GENERIC, [%'SSTM5079_144x91.nc'%], [%dnl
ifelse(FV_TYPE, DYCORE,  [%'SSTM5079_144x91.nc'%])dnl
%])dnl
%])dnl
%])dnl
%])dnl

define([%FV_SSTCYC%], [%dnl
ifelse(FV_TYPE, MONTHLY, [%.false.%], [%dnl
ifelse(FV_TYPE, WEEKLY,  [%.false.%], [%dnl
ifelse(FV_TYPE, GENERIC, [%.false.%], [%dnl
ifelse(FV_TYPE, DYCORE,  [%.false.%])dnl
%])dnl
%])dnl
%])dnl
%])dnl

define([%FV_OZNDATA%], [%dnl
ifelse(FV_TYPE, MONTHLY, [%'o3.amip2_uars_fub_91x55.nc'%], [%dnl
ifelse(FV_TYPE, WEEKLY,  [%'o3.amip2_uars_fub_181x55.nc'%],    [%dnl
ifelse(FV_TYPE, GENERIC, [%'o3.amip2_uars_fub_91x55.nc'%], [%dnl
ifelse(FV_TYPE, DYCORE,  [%'o3.amip2_uars_fub_91x55.nc'%])dnl
%])dnl
%])dnl
%])dnl
%])dnl

define([%FV_H2ODATA%], [%dnl
ifelse(FV_TYPE, MONTHLY, [%'RandelH2O_91x25.bin'%], [%dnl
ifelse(FV_TYPE, WEEKLY,  [%'RandelH2O_181x25.bin'%],    [%dnl
ifelse(FV_TYPE, GENERIC, [%'RandelH2O_91x25.bin'%], [%dnl
ifelse(FV_TYPE, DYCORE,  [%'RandelH2O_91x25.bin'%])dnl
%])dnl
%])dnl
%])dnl
%])dnl

define([%FV_SRFDATA%], [%dnl
ifelse(FV_TYPE, MONTHLY, [%'surf_r.data_144x91'%], [%dnl
ifelse(FV_TYPE, WEEKLY,  [%'surf_r.data_288x181'%],    [%dnl
ifelse(FV_TYPE, GENERIC, [%'surf_r.data_144x91'%], [%dnl
ifelse(FV_TYPE, DYCORE,  [%'surf_r.data_144x91'%])dnl
%])dnl
%])dnl
%])dnl
%])dnl

define([%FV_CCM3%], [%dnl
ifelse(FV_TYPE, MONTHLY, [%.true.%], [%dnl
ifelse(FV_TYPE, WEEKLY,  [%.true.%],    [%dnl
ifelse(FV_TYPE, GENERIC, [%.true.%],   [%dnl
ifelse(FV_TYPE, DYCORE, [%.false.%])dnl
%])dnl
%])dnl
%])dnl
%])dnl

define([%FV_NSPLIT%], [%dnl
ifelse(FV_TYPE, MONTHLY, [%4%], [%dnl
ifelse(FV_TYPE, WEEKLY,  [%8%],    [%dnl
ifelse(FV_TYPE, GENERIC, [%4%],   [%dnl
ifelse(FV_TYPE, DYCORE, [%4%])dnl
%])dnl
%])dnl
%])dnl
%])dnl


dnl --------------------
dnl  Main template body
dnl --------------------

define([%FV_TMPL%], [%dnl
#!/bin/csh -f

ifelse(FV_ARCH, SGI, [%dnl
# ------------------------------
changecom('/*', '*/')dnl
#PBS -l ncpus=NCPUGEN
changecom()dnl
#PBS -l walltime=12:00:00
#PBS -l mem=4gb
#PBS -S /bin/csh
#PBS -V
#PBS -j eo
# ------------------------------

%], [%dnl
ifelse(FV_ARCH, IBM, [%dnl
ifelse(FV_HOST, f01n, [%dnl
#
# NCCS IBM
#
#
#@ job_type = parallel
#@ job_name = fvcam
#@ output = $(job_name).out.$(jobid)
#@ error =  $(job_name).err.$(jobid)
#@ class = poe
#@ tasks_per_node = 1
changecom('/*', '*/')dnl
#@ node = FV_NODE
changecom()dnl
#@ wall_clock_limit = 0:30:00
#@ network.MPI = css0,shared,us
#@ queue
%], [%dnl
#
# LLNL IBM
#
#PSUB -tM 2:00
#PSUB -b climate
#PSUB -mb
#PSUB -me
#PSUB -e batch.SUBCASE.err
#PSUB -o batch.SUBCASE.out
#
# NERSC/ORNL/NCAR IBM
#
#@ job_name = SUBCASE
#@ job_type = parallel
#@ output = batch.SUBCASE.out
#@ error = batch.SUBCASE.err
#@ class = regular
#@ tasks_per_node = 1
#@ node = 1
#@ wall_clock_limit = 2:00:00
#@ notification    = always
#@ network.MPI = css0,shared,us
#@ environment     = COPY_ALL
#@ queue
%])dnl
%], [%dnl
ifelse(FV_ARCH, LINUX, [%dnl

%], [%dnl
ifelse(FV_ARCH, DEC, [%dnl
#PSUB -tM 2:00
#PSUB -b climate
#PSUB -mb
#PSUB -me
#PSUB -e batch.SUBCASE.err
#PSUB -o batch.SUBCASE.out

%], [%dnl
ifelse(FV_ARCH, CRAY, [%dnl
#QSUB -me
#QSUB -e batch.SUBCASE.err
#QSUB -o batch.SUBCASE.out
#QSUB -s /bin/csh

%], [%dnl
ifelse(FV_ARCH, CRAY_T3E, [%dnl
#QSUB -l mpp_p=16
#QSUB -l mpp_t=3600
#QSUB -l p_mpp_t=3600
#QSUB -me
#QSUB -e batch.SUBCASE.err
#QSUB -o batch.SUBCASE.out
#QSUB -s /bin/csh

%])dnl
%])dnl
%])dnl
%])dnl
%])dnl
%])dnl
#
# ... Is this INTERACTIVE or BATCH?
#
 if (! $?ENVIRONMENT) then
   setenv ENVIRONMENT INTERACTIVE
 endif
 echo "Environment variable defined as $ENVIRONMENT"
 if ("$ENVIRONMENT" == "INTERACTIVE") then
   set jobname = $[%%]0
 endif

ifelse(FV_ARCH, SGI, [%dnl
 if ("$ENVIRONMENT" == "BATCH") then
   set jobname = $PBS_JOBNAME
 endif

%], [%dnl
ifelse(FV_ARCH, IBM, [%dnl
 if ("$ENVIRONMENT" == "BATCH") then
   set jobname = dycore.j
 endif

%], [%dnl
ifelse(FV_ARCH, LINUX, [%dnl

%], [%dnl
ifelse(FV_ARCH, DEC, [%dnl
 if ("$ENVIRONMENT" == "BATCH") then
   set jobname = dycore.j
 endif

%], [%dnl
ifelse(FV_ARCH, CRAY, [%dnl
 if ("$ENVIRONMENT" == "BATCH") then
   set jobname = dycore.j
 endif

%], [%dnl
ifelse(FV_ARCH, CRAY_T3E, [%dnl
 if ("$ENVIRONMENT" == "BATCH") then
   set jobname = dycore.j
 endif

%])dnl
%])dnl
%])dnl
%])dnl
%])dnl
%])dnl

#
# -----------------------------------------------------
# Usage:
#        fvcam.j      :for normal restart run (default)
#        init.j       :for init run
#        branch.j     :for branch run
# -----------------------------------------------------
#
 set Restart = 1                # Normal restart run
 switch (${jobname})
 case init.j:
   set Restart = 0              # Initialization run
 breaksw
 case branch.j:
   set Restart = 2              # Branch run
 breaksw
 default:
 endsw
#
# ##############################
# #  Start User Configuration  #
# ##############################
#
# ... System Configuration
#
ifelse(FV_ARCH, SGI, [%dnl
 set archtype = "sgi"
#from CCM4
limit stacksize   unlimited                # KB
setenv OMP_DYNAMIC FALSE
setenv _DSM_PLACEMENT FIRST_TOUCH
setenv _DSM_WAIT SPIN
setenv MPC_GANG    OFF
setenv MP_SLAVE_STACKSIZE 40000000
setenv TRAP_FPE "UNDERFL=FLUSH_ZERO; OVERFL=ABORT,TRACE; DIVZERO=ABORT,TRACE"

%], [%dnl
ifelse(FV_ARCH, IBM, [%dnl
#setenv XLFRTEOPTS namelist=OLD
 setenv XLSMPOPTS "stack = 0x04000000"       # bytes
 setenv MP_LABELIO no
 set archtype = "ibm"

%], [%dnl
ifelse(FV_ARCH, LINUX, [%dnl
 limit stacksize   unlimited
# setenv MPSTKZ     32M
 set archtype = "linux"

%], [%dnl
ifelse(FV_ARCH, DEC, [%dnl
 setenv KMP_STACKSIZE 64000000             # bytes
 setenv MP_STACK_SIZE 64000000             # bytes
 set archtype = "dec"

%], [%dnl
ifelse(FV_ARCH, CRAY, [%dnl
 set archtype = "cray"

%], [%dnl
ifelse(FV_ARCH, CRAY_T3E, [%dnl
 set archtype = "cray"

%])dnl
%])dnl
%])dnl
%])dnl
%])dnl
%])dnl

 setenv MPI        FV_MPI
 setenv AGCM_N_PROCESSES      FV_N_MPI
 setenv AGCM_N_THREADS_PER_PROCESS      FV_N_SMP
 setenv BaseDir    FV_BASEDIR
 setenv VER        FV_VER
 setenv CASE       FV_CASE
 set CASEID = `echo $CASE | awk '{printf "%s",[%substr($%][%1,1,3)%]}' - `

dnl
dnl
ifelse(FV_TYPE, DYCORE,, [%dnl
dnl  if ($MPI != 0) then
dnl   echo ""
dnl   echo "MPI not supported for this configuration - disabling"
dnl   echo ""
dnl   setenv MPI 0
dnl   setenv AGCM_N_PROCESSES 1
dnl  endif
%]) dnl
dnl
dnl

#
 set EXE         = FV_EXE
 set LOG         = ${EXE}.log
 set MSWRITE     = mswrite          #for CCM4
 set MSREAD      = msread
 set FFC_RUN     = ""
 ifelse(FV_ARCH, LINUX, [%dnl
#   set FFC_RUN   = "-Wl,-T"
%]) dnl
#
 set SourceDir   = ${BaseDir}/${VER}
 set FVData      = FV_DATA
 set LMData      = FV_LMDATA
 set FVNCDATA    = FV_NCDATA
 set FVBNDTVS    = FV_BNDTVS
 set FVCLMSURF   = FV_CLM_SURF
 set s_name      = sstsice
# set s_res       = 288x181
 set s_freq      = weekly
#
 set cam_rest_pfile = "./cam.rpointer"
 set lsm_rest_pfile = "./lsm.rpointer"
 set lnd_rest_pfile = "./lnd.rpointer"

#
# Set Work directory according to architecture
#
ifelse(FV_ARCH, SGI, [%dnl
 set WorkDir0     = FV_SCRATCH
 set WorkDir1     = ${WorkDir0}
%], [%dnl
ifelse(FV_ARCH, IBM, [%dnl
ifelse(FV_MACH, seaborg, [%dnl
 set WorkDir0     = /scratch
%], [%dnl
 set WorkDir0     = FV_SCRATCH
%])dnl
ifelse(FV_HOST, f01n, [%dnl 
  set WorkDir0    = FV_SCRATCH
%])dnl
 set WorkDir1     = ${WorkDir0}/${LOGNAME}/${VER}/${CASE}
%], [%dnl
ifelse(FV_ARCH, LINUX, [%dnl
 set WorkDir0     = FV_SCRATCH
 set WorkDir1     = ${WorkDir0}
%], [%dnl
ifelse(FV_ARCH, DEC, [%dnl
 set WorkDir0     = FV_SCRATCH
 set WorkDir1     = ${WorkDir0}/${LOGNAME}/${VER}/${CASE}
%], [%dnl
ifelse(FV_ARCH, CRAY, [%dnl
 set WorkDir0     = /usr/tmp
 set WorkDir1     = ${WorkDir0}/${LOGNAME}/Work/${VER}/${CASE}
%], [%dnl
ifelse(FV_ARCH, CRAY_T3E, [%dnl
 set WorkDir0     = /usr/tmp
 set WorkDir1     = ${WorkDir0}/${LOGNAME}/Work/${VER}/${CASE}
%])dnl
%])dnl
%])dnl
%])dnl
%])dnl
%])dnl

 set WorkDir      = ${WorkDir1}/SUBCASE
 
#
# Postprocessing:
#
#     Interpolate diagnostics output from eta to pressure 
#     coordinates, with an option to save the interpolated
#     dataset (at pressure levels) instead of the original
#     dataset (at eta levels). SavePres=1 implies that the
#     eta dataset will not be archived to the MSS.
#
ifelse(FV_ARCH, SGI, [%dnl
 set PostProc    = 1     # 1=Yes,0=No
%],[%dnl
 set PostProc    = 0            # 1=Yes,0=No
%])dnl

 set SavePres    = 0            # 1=Yes,0=No

ifelse(FV_ARCH, SGI, [%dnl
 set SaveToSilo  = FV_SaveToSilo    # 1=Yes for a restart run
 set SaveToMSS   = FV_SaveToMSS     # 1=Yes,0=No
%],[%dnl
 set SaveToSilo  = FV_SaveToSilo
 set SaveToMSS   = FV_SaveToMSS     # 1=Yes,0=No
%])dnl

 set MSSHost     = helios1
ifelse(FV_STAGE, USESTAGE, [%dnl
#
 set UseStage    = 1            # 1=Yes,0=No
 set StageHost   = dixon0
 set StageDisk   = stage2
%])dnl
#
 set SpinUp      = 0            # 1=Yes,0=No
 set Benchmark   = 0            # 1=Yes,0=No
 set CleanSilo   = 0            # 1=Yes,0=No
 set RenHtapes   = 0            # 1=Yes,0=No



if ($SaveToSilo == 1 || $SaveToMSS == 1) then
   set    IRT = 1
   set   RIRT = 1
else
   set    IRT = 0
   set   RIRT = 0
endif



#
 set MaxCount    = 1
#
# ... Model Configuration
#
ifelse(FV_TYPE, MONTHLY, [%dnl
dnl ---------
dnl  monthly 
dnl ---------
 set Days   = (31 28 31 30 31 30 31 31 30 31 30 31)
 set Months = (Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec)
 
 set Year   = YEAR           # model year
 set Mon    = MONTH          # model month
 set it     = 3

 set NYMDE  = 19990101       # yyyymmdd, ending date

 set count  = 1
 while ( ${count} <= ${MaxCount} )       # Time Loop (monthly)

 set Year   = `echo ${Year} | awk '{printf "%4.4d", $[%%]1}' -`
 set Mon    = `echo ${Mon}  | awk '{printf "%2.2d", $[%%]1}' -`
 
 set NYMDE1 = "${Year}${Mon}00"
 echo "NYMDE1 $NYMDE1"
 if ( "${NYMDE1}" >= "${NYMDE}" ) then
   echo "Simulation already completed, stop."
   exit 0
 endif

    
#
# Check if $Year is a leap year 
#
  set Leap_Year = 0
  if ( "$Year" >= 1900 ) then
     if ( "`expr $Year \% 4`" == 0 ) then
        if ( "`expr $Year \% 100`" == 0 ) then
           if ( "`expr $Year \% 400`" == 0 ) then
               set Leap_Year = 1
           endif
        else
            set Leap_Year = 1
        endif
     endif
  endif
 if ( "$Leap_Year" == 1 ) then
#   set Days   = (31 29 31 30 31 30 31 31 30 31 30 31) #turned off 
   set Days   = (31 28 31 30 31 30 31 31 30 31 30 31)
 endif
#
 set mmsave = ${Mon}
 set days   = ${Days[${Mon}]}
 set days4  = `expr 4 \* ${days}`
 set s_date = y${Year}
 set s_res  = 144x91
%], [%dnl
ifelse(FV_TYPE, WEEKLY, [%dnl
dnl --------
dnl  weekly 
dnl --------
 set Days   = (31 28 31 30 31 30 31 31 30 31 30 31)
 set Months = (Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec)
 
 set Year   = YEAR           # model year
 set Mon    = MONTH          # model month
 set Week   = WEEK           # week of the month
 set it     = 3

 set count  = 1
 while ( ${count} <= ${MaxCount} )       # Time Loop (weekly)

 set Year   = `echo ${Year} | awk '{printf "%4.4d",$[%%]1}' -`
 set Mon    = `echo ${Mon}  | awk '{printf "%2.2d",$[%%]1}' -`
#
# Check if $Year is a leap year
#
  set Leap_Year = 0
  if ( "$Year" >= 1900 ) then
     if ( "`expr $Year \% 4`" == 0 ) then
        if ( "`expr $Year \% 100`" == 0 ) then
           if ( "`expr $Year \% 400`" == 0 ) then
               set Leap_Year = 1
           endif
        else
            set Leap_Year = 1
        endif
     endif
  endif
 if ( "$Leap_Year" == 1 ) then
#   set Days   = (31 29 31 30 31 30 31 31 30 31 30 31) #turned off
   set Days   = (31 28 31 30 31 30 31 31 30 31 30 31)
 endif
#
 set mmsave = ${Mon}
 set wksave = ${Week}
 set s_date = y${Year}
 set s_res  = 288x181


 if (${Week} == 1 ) then
   cp /dev/null ${SourceDir}/${CASE}/wklists.${Mon}.${Year}
 endif

 if (${Week} <= 4) then
   set days = 7
 else
   @ days   = ${Days[${Mon}]} - (${Week} - 1) * 7
 endif
%], [%dnl
ifelse(FV_TYPE, GENERIC, [%dnl
 set count  = 1
 while ( ${count} <= ${MaxCount} )       # Time Loop
%], [%dnl
ifelse(FV_TYPE, DYCORE, [%dnl
 set count  = 1
 while ( ${count} <= ${MaxCount} )       # Time Loop
%])dnl
%])dnl
%])dnl
%])dnl



ifelse(FV_ARCH, SGI, [%dnl
 set Workbase = `echo ${WorkDir0}|cut -d'/' -f2`
 set Versilo = ${VER}
%],[%dnl
 set Workbase = `echo ${WorkDir0}|cut -b 2-`
 set Versilo = Silo/${VER}
%])dnl

 if (! ${?Silo}) then
ifelse(FV_TYPE, MONTHLY, [%dnl
#08/23/200   set Silo   = "/${Workbase}/${LOGNAME}/${Versilo}/${CASE}/${Year}/${Months[${Mon}]}"
   set Silo   = "/${Workbase}/${LOGNAME}/${Versilo}/${CASE}"
   set SavDir = "/${Workbase}/${LOGNAME}/${VER}/${CASE}/${Year}"
%], [%dnl
ifelse(FV_TYPE, WEEKLY, [%dnl
   set Silo   = "/${Workbase}/${LOGNAME}/${Versilo}/${CASE}/${Year}/${Months[${Mon}]}/${Week}"
   set SavDir = "/${Workbase}/${LOGNAME}/${VER}/${CASE}/${Year}/${Months[${Mon}]}"
%], [%dnl
ifelse(FV_TYPE, [%GENERIC%], [%dnl
   set Silo   = "/${Workbase}/${LOGNAME}/${Versilo}/${CASE}"
   set SavDir = "/${Workbase}/${LOGNAME}/${VER}"
%], [%dnl
ifelse(FV_TYPE, [%DYCORE%], [%dnl
   set Silo   = "/${Workbase}/${LOGNAME}/${Versilo}/${CASE}"
   set SavDir = "/${Workbase}/${LOGNAME}/${VER}"
%])dnl
%])dnl
%])dnl
%])dnl
 endif


ifelse(FV_TYPE, MONTHLY, [%dnl
 set atm_parm = "ccm4.namelist.SUBCASE"
 if ( "$it" == 0 ) then
    set nsrest = 0
    set initmod = 'arbitrary initialization'
 else
    set nsrest = $it
    set initmod = ' '
 endif
dnl
dnl 05/30/20001
dnl
set NREVSN_HOME=""
set camr_year=""
set camr_file=""
set lsmr_file=""
set lndr_file=""
if ( -e ${SourceDir}/${CASE}/$cam_rest_pfile ) then
 set ffile=`head -1 ${SourceDir}/${CASE}/$cam_rest_pfile`
 set camr_file=`basename $ffile`
 #set camr_year=`echo $camr_file|awk '{printf "%s",[%substr($%][%1,6,4)%]}' - `
 set camr_year=`echo $camr_file| cut -d '.' -f2 | cut -d '-' -f1`
if ($SaveToSilo == 1 || $SaveToMSS == 1) then
 set NREVSN_HOME="${Silo}/rest/${camr_year}"
else
 set NREVSN_HOME="."
endif
 echo $NREVSN_HOME
endif
if ( -e ${SourceDir}/${CASE}/$lsm_rest_pfile ) then
 set ffile=`head -1 ${SourceDir}/${CASE}/$lsm_rest_pfile`
 set lsmr_file=`basename $ffile`
endif
if ( -e ${SourceDir}/${CASE}/$lnd_rest_pfile ) then
 set ffile=`head -1 ${SourceDir}/${CASE}/$lnd_rest_pfile`
 set lndr_file=`basename $ffile`
endif

%], [%dnl   else
 set atm_parm = "ccm4.namelist.SUBCASE"
  set nsrest = 0
  set initmod = 'arbitrary initialization'
%])dnl   End of ifelse Monthly

dnl Year not defined in dycore.j, set s_date=y${Year}
dnl
dnl create namelists
dnl
 cat >! ${SourceDir}/${CASE}/$atm_parm << EOF
ifelse(FV_DECOMP, twod, [%dnl
 &mprun2d
 npr_yz      = FV_NPR_YZ
 /
%])dnl  End of ifelse FV_DECOMP
 &CAMEXP
ifelse(FV_TYPE, MONTHLY, [%dnl
 caseid      = '${CASE}.${Year}-${Mon}'
%], [%dnl   else
 caseid      = '${CASEID}'
%])dnl   End of ifelse MONTHLY
 ctitle      = '${VER}_${CASE}'
 ncdata      = '${FVData}/inic/fv/${FVNCDATA}'
 bndtvs      = '${FVData}/sst/${FVBNDTVS}'
 bndtvo      = '${FVData}/ozone/noaao3.1990.21999.nc'
 bndtvg      = '${FVData}/ggas/noaamisc.r8.nc'
 absems_data ='${FVData}/rad/abs_ems_factors_fastvx.052001.nc'
 rest_pfile  = '$cam_rest_pfile'
 archive_dir = '$Silo'
 iradsw      =  -1    
 iradlw      =  -1    
 iradae      =  -12    
 dtime       =  FV_DTIME
 use_eta     = .true.
 mss_irt     =  $IRT
 ndens       = 2,2,2
ifelse(FV_TYPE, MONTHLY, [%dnl
 mfilt       = 1,${days},${days4}
 nelapse     = -${days}
 nhtfrq      = 0,-24,-6
 NREVSN      = '${NREVSN_HOME}/${camr_file}'
%], [%dnl   else
 mfilt       = 1,1
 nelapse     = -1
 nhtfrq      =  0
%])dnl   End of ifelse MONTHLY
 nsrest      = ${nsrest}
 iyear_ad    = 1950
 readtrace   = .false.
 trace_gas   = FV_TRACE_GAS
 adiabatic   = .false.
 ideal_phys  = .false.
 inithist    = 'MONTHLY'
 fincl3      = 'V:I', 'U:I','T:I','Q:I','CWAT:I','OMEGA:I','PS:I' 
 nsplit      = 0
 nrefrq      = 1
 print_step_cost = .F.
 reset_csim_iceprops = .T. 
 /
 &lsmexp
 archive_dir = '$Silo'
 mksrfdir    = '$FVData'
 lsmgeo      = FV_LSMGEO
 finidat     = '$initmod'
 rest_pfile  = '$lsm_rest_pfile'
ifelse(FV_TYPE, MONTHLY, [%dnl
 NREVSN      = '${NREVSN_HOME}/${lsmr_file}'
%])dnl   End of ifelse MONTHLY
 /
 &clmexp
 archive_dir    = '$Silo'
 fsurdat        = '$LMData/srfdata/cam/${FVCLMSURF}'
 fpftcon        = '$LMData/pftdata/pft-physiology'
 mksrf_fvegtyp  = '$LMData/rawdata/mksrf_pft.nc'
 mksrf_fsoitex  = '$LMData/rawdata/mksrf_soitex.10level.nc'
 mksrf_fsoicol  = '$LMData/rawdata/mksrf_soicol_clm2.nc'
 mksrf_flanwat  = '$LMData/rawdata/mksrf_lanwat.nc'
 mksrf_furban   = '$LMData/rawdata/mksrf_urban.nc'
 mksrf_fglacier = '$LMData/rawdata/mksrf_glacier.nc'
 mksrf_flai     = '$LMData/rawdata/mksrf_lai.nc'
ifelse(FV_TYPE, MONTHLY, [%dnl
 nrevsn         = '${NREVSN_HOME}/${lndr_file}'
%])dnl   End of ifelse MONTHLY
 /
EOF


#
# ############################
# #  End User Configuration  #
# ############################
#
# ... Configurations derived from user section
#

ifelse(FV_ARCH, SGI, [%dnl
 setenv MP_SET_NUMTHREADS   ${AGCM_N_THREADS_PER_PROCESS}
 setenv OMP_NUM_THREADS     ${AGCM_N_THREADS_PER_PROCESS}
 set stacksize = `limit stacksize | awk '{print $[%%]2}' -`
 if (${stacksize} == "unlimited" || ${stacksize} > 4000000) then
##error in running ccm4   setenv MP_STACK_OVERFLOW OFF
    echo " "
 endif

%], [%dnl
ifelse(FV_ARCH, IBM, [%dnl
 setenv XLSMPOPTS "$XLSMPOPTS :  parthds=${AGCM_N_THREADS_PER_PROCESS} : schedule=guided"
 if ("$ENVIRONMENT" == "BATCH") then
ifelse(FV_MACH, seaborg, [%dnl
   setenv MP_RMPOOL 2
%], [%dnl
   setenv MP_RMPOOL 1
%])dnl
 else
ifelse(FV_MACH, seaborg, [%dnl
   setenv MP_RMPOOL 2
%], [%dnl
   setenv MP_RMPOOL 0
%])dnl
   setenv MP_PROCS ${AGCM_N_PROCESSES}
   if ( "$AGCM_N_PROCESSES" == 0 ) then
      setenv MP_PROCS 1
   endif
 endif

%], [%dnl
ifelse(FV_ARCH, LINUX, [%dnl
 setenv OMP_NUM_THREADS     ${AGCM_N_THREADS_PER_PROCESS}
 limit stacksize unlimited
 set stacksize = `limit stacksize | awk '{print $[%%]2}' -`
 if (${stacksize} == "unlimited" || ${stacksize} > 4000000) then
##error in running ccm4   setenv MP_STACK_OVERFLOW OFF
    echo " "
 echo "setenv OMP_NUM_THREADS     ${AGCM_N_THREADS_PER_PROCESS}"  > ~/.cshrc4mpi
 echo "limit stacksize unlimited          " >>  ~/.cshrc4mpi
 endif

%], [%dnl
ifelse(FV_ARCH, DEC, [%dnl
 setenv OMP_SCHEDULE GUIDED
 setenv OMP_NUM_THREADS ${AGCM_N_THREADS_PER_PROCESS}

%], [%dnl
ifelse(FV_ARCH, CRAY, [%dnl

%], [%dnl
ifelse(FV_ARCH, CRAY_T3E, [%dnl

%])dnl
%])dnl
%])dnl
%])dnl
%])dnl
%])dnl

 echo "Silo    --> ${Silo}"
 if ("${Silo}" == "${SourceDir}/${CASE}") then
   set Samesilo = 1
   echo "Silo is coincident with Case area of Source Directory"
 else
   set Samesilo = 0
   if (${archtype} == "cray") then
#     rm -fr ${Silo}
      echo 'hi'
   endif
   mkdir -p ${Silo}
 endif

#
#


 if (${SaveToMSS}) then
ifelse(FV_STAGE, USESTAGE, [%dnl
   set LocalStage = 1
%])dnl
   if (! ${?MSSPath}) then
     set MSSPath = ${VER}/${CASE}
dnl  ifelse(FV_TYPE, MONTHLY, [%/${Year}%], [%dnl
dnl  ifelse(FV_TYPE, WEEKLY,  [%/${Year}/${Months[${Mon}]}%], [%dnl
dnl  ifelse(FV_TYPE, GENERIC, [%%], [%dnl
dnl  ifelse(FV_TYPE, DYCORE, [%%])dnl
dnl  %])dnl
dnl  %])dnl
dnl  %])
   endif
# --- Only add -hp suffix if it is not already there
   if (! `echo ${MSSHost} | grep '\-hp$' | wc -c`) then
     set MSSHost = ${MSSHost}-hp
   endif
   set MSS = ${MSSHost}:${MSSPath}
   echo "MSS Base     --> ${MSS}"
 else
ifelse(FV_STAGE, USESTAGE, [%dnl
   set UseStage  = 0
%])dnl
   set MSSPath   =
   set MSS       =
 endif

ifelse(FV_STAGE, USESTAGE, [%dnl
 if (${UseStage}) then
   if (! ${?StagePath}) then
     set StagePath = /${StageDisk}/${LOGNAME}/${VER}/${CASE}dnl
ifelse(FV_TYPE, MONTHLY, [%/${Year}/${Months[${Mon}]}%], [%dnl
ifelse(FV_TYPE, WEEKLY,  [%/${Year}/${Months[${Mon}]}/${Week}%], [%dnl
ifelse(FV_TYPE, GENERIC, [%%], [%dnl
ifelse(FV_TYPE, GENERIC, [%%])dnl
%])dnl
%])dnl
%])
   endif
   set Stage = ${StageHost}-hp:${StagePath}
   echo "Stage   --> ${Stage}"
   if ("${StageHost}" == "`/usr/bsd/hostname`") then
     set LocalStage = 1
     set StagePath  = ${Silo}
   else
     set LocalStage = 0
   endif
 else
   set UseStage = 0
 endif
%])dnl
#
# ... Set up working directory and copy/link all the restart files along
#     with model executable, utilities, namelist and diagnostic table
#     file to this directory (if necessary).
#     The simulation is carried out here in the working directory.
#
#     See if distinct working directory is called for
#
 if ("${WorkDir}" == "${SourceDir}/${CASE}") then
   set Samedir = 1
   echo "Working in Case area of Source Directory"
 else
   set Samedir = 0
   if (${archtype} == "cray") then
     rm -fr ${WorkDir}
   endif
   mkdir -p ${WorkDir}
 endif
#
# Use of non-distinct working directory not yet fully supported (AAM)
#
 echo "WorkDir --> ${WorkDir}"
 echo " "
 cd ${WorkDir}
##  echo "WorkDir --> ${Silo} in Silo"
##  echo " "
##  cd ${Silo}


#
#
#


cat >! ${SourceDir}/${CASE}/$MSWRITE << EOF
#!/bin/csh -f
set IO = \$[%%]0
set idx=\$[%%]#argv
set pathname=\$argv[\$idx]
set filename=\`basename \$pathname \`
set tmp = \`dirname \$pathname \`
set dir = \`basename \$tmp \`          #hist or rest
#
set h2 = \` echo \$filename | cut -d'.' -f2 - \`
set h3 = \` echo \$filename | cut -d'.' -f3 - \`
set h4 = \` echo \$filename | cut -d'.' -f4 - \`
#
set f4 = \` echo \$filename | awk '{printf "%s",[%substr(\$%][%1,1,4)%]}' - \`
switch (\$f4)
  case cami:
  case lsmi:  # not good for lsmi_arbini.nc
      if ("\$filename" == "lsmi_arbini.nc") then
       set dir = "\${dir}/lsm"
      else
       set dir = "\${dir}/init"
       set yr=\`echo \$filename|awk -F'_' '{printf "%s", [%substr(\$%][%2,1,4)%]}'\`
       set dir = "\${dir}/\${yr}"
      endif
       breaksw
  case lsm_:
       set dir = "\${dir}/lsm"
       breaksw
  case clmr:             #added 10/12/2001
  case lsmh:             #hist
  case camr:             #rest
  case lsmr:
      if ("\$filename" == "lsmh_timcon.nc") then
       set dir = "\${dir}/lsm"
      else
       set yr=\`echo \$filename|awk -F'_' '{printf "%s", [%substr(\$%][%2,1,4)%]}'\`
       set dir = "\${dir}/\${yr}"
      endif
       breaksw
  case "ha*":
       set dir = "\${dir}/ha"
       breaksw
  case "hb*":
       set dir = "\${dir}/hb"
       breaksw
  case "hc*":
       set dir = "\${dir}/hc"
       breaksw
  case [0-9][0-9][0-9][0-9]:
       set yr="\$f4"
       set dir = "\${dir}/\${yr}"
       breaksw
  default:
       #echo "dir not changed"
       breaksw
endsw
#
set yr = \` echo \$h2 | awk '{printf "%s",[%substr(\$%][%1,1,4)%]}' - \` 
set dir = "\${dir}/\${yr}"
#
set Silo = $Silo
set SaveToMSS  = $SaveToMSS
set SaveToSilo = $SaveToSilo
set MSS        = $MSS
set MSSHost    = $MSSHost
set MSSPath    = $MSSPath
echo "To: \$Silo/\${dir}"
set IO = \`basename \$IO\`
if ( "\$IO" == "mswrite" ) then
   if (\$SaveToSilo == 1 ) then
      echo "mswrite \$filename to \$Silo/\${dir} in c-shell"
      mkdir -p \$Silo/\${dir}
      cp \$filename               \$Silo/\${dir}
   endif
   if (\$SaveToMSS == 1 ) then
      echo "mswrite \$filename to \${MSS}/\${dir} in c-shell"
      ssh \${MSSHost}  "mkdir -p \${MSSPath}/\${dir}"
#      scp -v \$filename            \${MSS}/\${dir}
      scp  \$filename            \${MSS}/\${dir}
   endif
else if ( "\$IO" == "msread" ) then
   echo "msread \$filename in c-shell"
   cp \$Silo/\${dir}/\$filename .
endif

exit 0

#scp to mass storage
EOF
chmod 755 ${SourceDir}/${CASE}/$MSWRITE

#
#
#
cp ${SourceDir}/${CASE}/$MSWRITE  ${SourceDir}/${CASE}/$MSREAD 


#
# copy or link, depending on architecture
#
ifelse(FV_ARCH, SGI, [%dnl
 /bin/cp -f ${SourceDir}/${CASE}/$EXE                     .
 /bin/cp -f ${SourceDir}/${CASE}/$atm_parm                .
 /bin/cp -f ${SourceDir}/${CASE}/$MSWRITE                 .
 /bin/cp -f ${SourceDir}/${CASE}/$MSREAD                  .
%],[%dnl
 ln -s -f ${SourceDir}/${CASE}/$EXE                     $EXE
 ln -s -f ${SourceDir}/${CASE}/$atm_parm                  .
 /bin/cp -f ${SourceDir}/${CASE}/$MSWRITE                 .
 /bin/cp -f ${SourceDir}/${CASE}/$MSREAD                  .
%])dnl

if  ( "$nsrest" == "1" || "$nsrest" == "3") then
  /bin/cp -f ${SourceDir}/${CASE}/${cam_rest_pfile}      . >& /dev/null
  if ( -e ${SourceDir}/${CASE}/${lsm_rest_pfile} ) then
  /bin/cp -f ${SourceDir}/${CASE}/${lsm_rest_pfile}      . >& /dev/null
  endif

#10/12/2001
ifelse(FV_TYPE, MONTHLY, [%dnl
  if (-e ~/lnd.${CASE}.${Year}-${Mon}.rpointer ) then
   cp -f ~/lnd.${CASE}.${Year}-${Mon}.rpointer ./$lnd_rest_pfile  >& /dev/null
  endif
%], [%dnl   else
  if (-e  ~/lnd.${CASE}.rpointer ) then
   cp -f ~/lnd.${CASE}.rpointer                ./$lnd_rest_pfile  >& /dev/null
  endif
%])dnl

  if ( -e ./$cam_rest_pfile ) then 
  set tmp=`head ./$cam_rest_pfile`
  set camr_file=`basename $tmp`
  set yr=`echo $camr_file|awk -F'_' '{printf "%s", [%substr($%][%2,1,4%])}'`
  set tmp=`dirname $tmp`
  set dir=`basename $tmp`
  set camr_file="${Silo}/${dir}/${yr}/${camr_file}"  
  echo "CAMR is $camr_file"
  endif

  set lsmr_file=""
  if ( -e ./$lsm_rest_pfile ) then
  set tmp=`head ./$lsm_rest_pfile`
  set lsmr_file=`basename $tmp`
  set yr=`echo $lsmr_file|awk -F'_' '{printf "%s", [%substr($%][%2,1,4%])}'`
  set tmp=`dirname $tmp`
  set dir=`basename $tmp`
  set lsmr_file="${Silo}/${dir}/${yr}/${lsmr_file}"
  echo "LSMR is $lsmr_file" 
  endif

  set clmr_file=""
  if ( -e ./$lnd_rest_pfile ) then
  set tmp=`head ./$lnd_rest_pfile`
  set clmr_file=`basename $tmp`
  set yr=`echo $clmr_file|awk -F'_' '{printf "%s", [%substr($%][%2,1,4%])}'`
  set tmp=`dirname $tmp`
  set dir=`basename $tmp`
  set clmr_file="${Silo}/${dir}/${yr}/${clmr_file}"
  echo "CLMR is $clmr_file"
  endif

 
  /bin/cp -f ${camr_file}                           .   >& /dev/null

  if ( "${lsmr_file}" != "" ) then
   if ( -e ${lsmr_file} ) then
    /bin/cp -f ${lsmr_file}                           .   >& /dev/null
   endif
  endif

  if ( "${clmr_file}" != "" ) then
   if ( -e ${clmr_file} ) then
    /bin/cp -f ${clmr_file}                           .   >& /dev/null
   endif
  endif
  /bin/cp -f ${camr_file}.A                         .   >& /dev/null
# /bin/cp -f ${lsmr_file}.A                         .
endif


ifelse(FV_FAKE, TRUE, [%dnl
   /bin/cp -f ${SourceDir}/${CASE}/v_fvccm3.sh             v_fvccm3.sh 
   /bin/cp -f ${SourceDir}/${CASE}/rst_stamp              rst_stamp
%]) dnl
 
#


ifelse(FV_FAKE, TRUE, [%dnl
  (time ./v_fvccm3.sh) |& tee $LOG
%]) dnl

set RUN_OPT_BEGIN =
set RUN_OPT_END =

ifelse(FV_ARCH, SGI, [%dnl
if ($MPI != 0) then
  set RUN_OPT_BEGIN = "mpirun -np $AGCM_N_PROCESSES"
endif

set RUN_OPT_BEGIN = "env MPC_GANG=off ${RUN_OPT_BEGIN}"

%], [%dnl
ifelse(FV_ARCH, IBM, [%dnl
set RUN_OPT_BEGIN = "poe"

%], [%dnl
ifelse(FV_ARCH, LINUX, [%dnl
  set RUN_OPT_END = $FFC_RUN
if ($MPI != 0) then
  set RUN_OPT_BEGIN = "mpirun -np $AGCM_N_PROCESSES"
endif


%], [%dnl
ifelse(FV_ARCH, DEC, [%dnl
if ("$ENVIRONMENT" == "BATCH") then
  set jobq = pbatch
else
  set jobq = pdebug
endif
if ($MPI != 0) then
  set RUN_OPT_BEGIN = "dmpirun -np $AGCM_N_PROCESSES"
# set RUN_OPT_BEGIN = "prun -p $jobq -N $AGCM_N_PROCESSES -n $AGCM_N_PROCESSES -c $AGCM_N_THREADS_PER_PROCESS"
endif

%], [%dnl
ifelse(FV_ARCH, CRAY, [%dnl

%], [%dnl
ifelse(FV_ARCH, CRAY_T3E, [%dnl
if ($MPI != 0) then
  set RUN_OPT_BEGIN = "mpprun -n $AGCM_N_PROCESSES"
endif

%])dnl
%])dnl
%])dnl
%])dnl
%])dnl
%])dnl

#
setenv PATH ${WorkDir}:${Silo}:${PATH}
#

( time ${RUN_OPT_BEGIN} ./$EXE ${RUN_OPT_END} < $atm_parm ) |&  tee $LOG

 set failed = $status

cp timing* ${SourceDir}/${CASE} >& /dev/null  #timing results


ifelse(FV_TYPE, MONTHLY, [%dnl
if ( "${RenHtapes}" == 1 ) then
   mv ${Silo}/hist/ha/ha0001.nc  ${Silo}/hist/ha/ha_${Year}-${Mon}.nc
   mv ${Silo}/hist/hb/hb0001.nc  ${Silo}/hist/hb/hb_${Year}-${Mon}.nc
endif
%])dnl   End of ifelse MONTHLY        

 if ( "${failed}" != 0 ) then
   echo "This simulation was failed with exit status = ${failed}, stop."
   exit (${failed})
 endif
#
#
#


# set tmp=`head ./$cam_rest_pfile`
# set cam_tmp=`basename $tmp`
# set tmp=`head ./$lsm_rest_pfile`
# set lsm_tmp=`basename $tmp`
# cat ./$cam_rest_pfile
# cat ./$lsm_rest_pfile
# echo "cam_tmp $cam_tmp lsm_tmp $lsm_tmp"

ifelse(FV_TYPE, MONTHLY, [%dnl
### cp $cam_rest_pfile ${Silo}/${cam_rest_pfile}_${Year}_${Mon}
### cp $lsm_rest_pfile ${Silo}/${lsm_rest_pfile}_${Year}_${Mon}
cp $cam_rest_pfile ${Silo}/${cam_rest_pfile} >& /dev/null
cp $lsm_rest_pfile ${Silo}/${lsm_rest_pfile} >& /dev/null
%], [%dnl   else
dnl
%])dnl   End of ifelse MONTHLY
# echo "./$cam_tmp" >$cam_rest_pfile
# echo "./$lsm_tmp" >$lsm_rest_pfile

cp -f $cam_rest_pfile  ${SourceDir}/${CASE} >& /dev/null
cp -f $lsm_rest_pfile  ${SourceDir}/${CASE} >& /dev/null

ifelse(FV_TYPE, MONTHLY, [%dnl
cp -f ~/lnd.${CASE}.${Year}-${Mon}.rpointer ${SourceDir}/${CASE}/$lnd_rest_pfile
%], [%dnl   else
cp -f ~/lnd.${CASE}.rpointer ${SourceDir}/${CASE}/$lnd_rest_pfile
dnl
%])dnl   End of ifelse MONTHLY


#
#

ifelse(FV_TYPE, MONTHLY, [%dnl
mkdir -p ${Silo}/namelists
mv $atm_parm ${Silo}/namelists/${atm_parm}_${Year}_${Mon}
%], [%dnl   else
mv $atm_parm ${Silo}/${atm_parm}
%])dnl   End of ifelse MONTHLY


#
#


cd ${SourceDir}/${CASE}

#
#
#

ifelse(FV_TYPE, MONTHLY, [%dnl
dnl
dnl --- If it's monthly job script 
dnl
       @ Mon ++
       @ Mon = ($Mon - 1) % 12 + 1
       if ($Mon == 1) @ Year ++
       set Mon  = `echo $Mon   | awk '{printf "%2.2d", $[%%]1}' -`
       set Year = `echo $Year  | awk '{printf "%4.4d", $[%%]1}' -`
       set month = `echo month | tr 'a-z' 'A-Z'`
       set year  = `echo year  | tr 'a-z' 'A-Z'`
       if (${count} == ${MaxCount}) then
         echo "Continue to re-sumbit next job."
         echo " "
         sed -e "s%set Mon    = ${month}%set Mon    = ${Mon}%1" \
             -e "s%set Year   = ${year}%set Year   = ${Year}%1" \
                monthly.j.tmpl >! ${CASE}_m.${Mon}
         chmod ug+x ${CASE}_m.${Mon}
         qsub ${CASE}_m.${Mon}
       endif
%], [%dnl
ifelse(FV_TYPE, WEEKLY, [%dnl
dnl
dnl --- If it's weekly job script 
dnl
       @ rest   = ${Days[${Mon}]} - ${Week} * 7
       if (${rest} > 0) then
         @ Week ++
       else
         set Week = 1                  # Next month
         @ Mon ++
         @ Mon = ($Mon - 1) % 12 + 1
         if ($Mon == 1) @ Year ++
       endif
       set Year  = `echo $Year  | awk '{printf "%4.4d", $[%%]1}' -`
       set Mon   = `echo $Mon   | awk '{printf "%2.2d", $[%%]1}' -`
       set year  = `echo year  | tr 'a-z' 'A-Z'`
       set month = `echo month | tr 'a-z' 'A-Z'`
       set week  = `echo week  | tr 'a-z' 'A-Z'`
       if (${count} == ${MaxCount}) then
         echo "Continue to re-sumbit next job."
         echo " "
         sed -e "s%set Mon    = ${month}%set Mon    = ${Mon}%1" \
             -e "s%set Year   = ${year}%set Year   = ${Year}%1" \
             -e "s%set Week   = ${week}%set Week   = ${Week}%1" \
                weekly.j.tmpl >! ${CASE}_w.${Mon}-${Week}
         chmod ug+x ${CASE}_w.${Mon}-${Week}
         qsub ${CASE}_w.${Mon}-${Week}
       endif
%], [%dnl
ifelse(FV_TYPE, GENERIC, [%dnl
dnl
dnl --- If it's generic jobs script fvcam.j
dnl
       if (${count} == ${MaxCount}) then
         echo "Continue to re-sumbit next job."
         echo " "
#         qsub ${jobname}
       endif
%], [%dnl
ifelse(FV_TYPE, DYCORE, [%dnl
dnl
dnl --- If it's dycore jobs script dycore.j
dnl
       if (${count} == ${MaxCount}) then
         echo "Continue to re-sumbit next job."
         echo " "
         qsub ${jobname}
       endif
%])dnl
%])dnl
%])dnl
%])dnl

ifelse(FV_STAGE, USESTAGE, [%dnl
# Set up directories in mass storage system

 if (${SaveToMSS}) then

   if (${UseStage}) then
     set RemoteHost = ${StageHost}-hp
     set RemotePath = ${StagePath}
   else
     set RemoteHost = ${MSSHost}
     set RemotePath = ${MSSPath}
   endif
   set Destination = ${RemoteHost}:${RemotePath}

   if (! ${LocalStage}) then
     echo  " --> Timing for staging ... "
     date 

#  Copy the following to mass storage system (or stage storage)

     ssh ${RemoteHost}  "mkdir -p ${RemotePath}" >& /dev/null
#     scp -v ${d_rst} ${p_rst} ${rout} ${Destination}
     if (${status}) then
       echo "scp failed, check the following files on ${Destination}"
       foreach f (${d_rst} ${p_rst} ${rout})
         echo " --> ${f}"
         end
       set CleanSilo = 0
     else
       if (${CleanSilo}) then
         /bin/rm -f ${d_rst} ${p_rst} ${rout}
       endif
     endif


     date

# clean-up

     if (${CleanSilo}) then
       cd ..
#       /bin/rm -fr ${Silo}
     else
       mail ${LOGNAME} << EOF
 *** This is a message sent from ${jobname}:

 Please remember to clean up ${Silo} on ${HOST}.
EOF
     endif   # CleanSilo

   endif     # LocalStage

 endif       # SaveToMSS
%])dnl



 
 @ count ++
 end          # End of Time Loop

 exit 0
%])dnl
dnl
divert(0)dnl
dnl
dnl ---------------------------------------------------
dnl
dnl -----------------------------
dnl  Include the definition file
dnl -----------------------------
include([%define.m4%])dnl
FV_TMPL[%%]dnl
dnl
dnl ---------------------------------------------------
dnl
