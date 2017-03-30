#!/bin/sh 
#
main(){
export TESTS AGCM_N_THREADS_PER_PROCESS AGCM_N_PROCESSES OMP USE_MPI1 USE_MPI2 USE_MLP CASE CASES Case SaveToSilo Scratch Ver
export COMM TMP Home LOG EXELOG JOBID Host WorkDir WorkDir1 Interactive
export SaveToMSS RunScript USER_FC  TESTDATA  IDS Sdatep1 Edate
TESTDATA=fvgcm-1_3_36
#CASE=z18
SaveToSilo=""
SaveToMSS=""
Scratch=""
Home=`pwd`
TMP=${Home}/tmp.log
IDS=${Home}/.ids
>$IDS
RunScript=${Home}/fvtest.j
LOG=log            #don't give full path
EXELOG=$LOG
CPRNC=$HOME/bin/cprnc
JOBID=""
USER_FC=""
Interactive=0
Host=`hostname |cut -d"." -f1 `
#
#VER=fvgcm     $bug
tmpver=`(cd ..; pwd)`
VER=`basename ${tmpver}`
Tag="${Home}/../CVS/Tag"
if [ -s ${Tag}  ]; then
   tmp=`grep ccm $Tag`
#   VER=`echo $tmp | sed 's/Nccm3/fvccm4/g'  `  #Bug
fi
#
>${LOG}.2
Sdatep1="00470902"
Edate="00470903"
case $Host in
  dycore)
    create_qsub
    if [ "$GCase" = "" ]; then
       CASE=a18
    else
       CASE=$GCase
    fi
    Sdatep1="20170124"
    Edate="20170125"
    TESTS="L1 L2 "
    Scratch='~/scratch/$VER/$CASE'
    Interactive=1
    SaveToSilo=0
    SaveToMSS=0
    USER_FC=lf95
    ;;
  kalnay|daley|mintz)
    TESTS="K1 K2 K3 "
    if [ "$GCase" = "" ]; then
       CASE=b55
    else
       CASE=$GCase
    fi
    if [ "$CASE" = b55 -o "$CASE" = c55 ]; then
       if [ "$GPhy" = 2 ]; then
          TESTS="K1 K2 K3 M1 M2 "     #lsm
       else
          TESTS="K1 K2 K3 "
       fi
    else
       TESTS="L1 L2 L3"
    fi
#    Scratch='/scratch1/$LOGNAME/$VER/$CASE'
#not defined please    EXELOG=""
    ;;
  tropic)
    if [ "$GCase" = "" ]; then
       CASE=b55
    else
       CASE=$GCase
    fi
    Scratch='~/scratch/$VER/$CASE'
    TESTS="T1 T2 T3"
    SaveToSilo=0
    SaveToMSS=0
#   Interactive=1
    ;;
  f01*)
    TESTS="F1 F2 F3"
    ;;
  *)
    TESTS=""
    ;;
esac
#
get_headers4RunScript
#
echo "=========================================="
count=1
for i in $TESTS
do
 OMP=0
 USE_MPI1=0
 USE_MPI2=0
 USE_MLP=0
 N_MLP=0
 AGCM_N_PROCESSES=0
 AGCM_N_THREADS_PER_PROCESS=0
 case $i in
   L1)
     AGCM_N_THREADS_PER_PROCESS=0
     AGCM_N_PROCESSES=1
     USE_MPI1=1
     ;;
   L2)
     AGCM_N_THREADS_PER_PROCESS=0
     AGCM_N_PROCESSES=2
     USE_MPI1=1
     ;;
   L3)
     AGCM_N_THREADS_PER_PROCESS=2
     AGCM_N_PROCESSES=0
     USE_MPI1=0
     ;;
   K1)
     AGCM_N_THREADS_PER_PROCESS=16
     AGCM_N_PROCESSES=0
     USE_MPI2=0
     ;;
   K2)
     AGCM_N_THREADS_PER_PROCESS=0
     AGCM_N_PROCESSES=16
     USE_MPI2=1
     ;;
   K3)
     AGCM_N_THREADS_PER_PROCESS=4
     AGCM_N_PROCESSES=4
     USE_MPI2=1
     ;;
   M1)
     AGCM_N_THREADS_PER_PROCESS=0
     N_MLP=16
     ;;
   M2)
     AGCM_N_THREADS_PER_PROCESS=4
     N_MLP=4
     ;;
   T1)
     AGCM_N_THREADS_PER_PROCESS=7
     AGCM_N_PROCESSES=0
     ;;
   T2)
     AGCM_N_THREADS_PER_PROCESS=0
     AGCM_N_PROCESSES=3
     ;;
   T3)
     AGCM_N_THREADS_PER_PROCESS=1
     AGCM_N_PROCESSES=4
     ;;
   F1)
     AGCM_N_THREADS_PER_PROCESS=4
     AGCM_N_PROCESSES=0
     ;;
   F2)
     AGCM_N_THREADS_PER_PROCESS=0
     AGCM_N_PROCESSES=3
     ;;
   F3)
     AGCM_N_THREADS_PER_PROCESS=1
     AGCM_N_PROCESSES=4
     ;;
   *)
     continue
     ;;
 esac

 if [ "$GPhy" = 2 ]; then
    str="lsm"
 else
    str="clm"
 fi

 Case=${CASE}${str}
 Ver="$VER"  #initialized for each run

 if [ "$AGCM_N_THREADS_PER_PROCESS" -gt 0 ]; then
     OMP=1
#    Ver="${Ver}_OMP${AGCM_N_THREADS_PER_PROCESS}"
     Case="${Case}_O${AGCM_N_THREADS_PER_PROCESS}"
 fi
 if [ "$AGCM_N_PROCESSES" -gt 0 ]; then
###    USE_MPI2=1
#    Ver="${Ver}_MPI${AGCM_N_PROCESSES}"
    Case="${Case}_M${AGCM_N_PROCESSES}"
 fi
 if [ "$N_MLP" -gt 0 ]; then
    USE_MLP=1
    AGCM_N_PROCESSES=${N_MLP}
#    Ver="${Ver}_MLP${AGCM_N_PROCESSES}"
    Case="${Case}_ML${AGCM_N_PROCESSES}"
 fi
#    echo "Configuring (CASE=$CASE, AGCM_N_THREADS_PER_PROCESS=$AGCM_N_THREADS_PER_PROCESS, AGCM_N_PROCESSES=$AGCM_N_PROCESSES)"
    echo "Configuring $Case "
#
    cd $Home
    run_configure
    mv ../conf.log conf.log-${Case}
#
    CASES="$CASES $Case"
#

    #echo "exit #===========Skip testing...===================" >>$RunScript


done    #end of "for" TESTS?
#
run_fvtest
#
#

}


run_configure(){

rm -f conf.log > $TMP 2>&1
                     #v------------/ not \
cd ${Home}/..
rm -f conf.log
# a bug
#
if [ "$AGCM_N_PROCESSES" = 0 ]; then
   AGCM_N_PROCESSES=1
fi
if [ "$AGCM_N_THREADS_PER_PROCESS" = 0 ]; then
   AGCM_N_THREADS_PER_PROCESS=1
fi
echo "AGCM_N_THREADS_PER_PROCESS=$AGCM_N_THREADS_PER_PROCESS \n /      
      AGCM_N_PROCESSES=$AGCM_N_PROCESSES \n /
      OMP=$OMP \n     /
      USE_MPI1=$USE_MPI1 \n     /
      USE_MPI2=$USE_MPI2 \n     /
      USE_MLP=$USE_MLP \n     /
      SaveToSilo=$SaveToSilo \n /
      SaveToMSS=$SaveToMSS \n /
      Scratch="$Scratch" \n /
      Ver=$Ver \n /
      Case=${Case}  \n/
      READRST=0     \n/
      Phy=${GPhy}   \n      /
      USER_FC=$USER_FC \n" | configure_fv.pl -f - \
      | awk '/Quick/,/fvgcm/ {print $0} ' > $TMP
#      | grep -v '='                     > $TMP
cd ${Home}
#empty string should be in the last
}

get_headers4RunScript()
{
echo "#!/bin/csh -f" >$RunScript 
case $Host in
   dycore)
cat >>$RunScript <<EOF
alias qsub ${Home}/qsub
alias qrls echo
EOF
     ;;
   kalnay|daley|mintz)
cat >>$RunScript <<EOF
# ------------------------------
#PBS -l ncpus=16
#PBS -l walltime=12:00:00
#PBS -l mem=4gb
#PBS -S /bin/csh
#PBS -V
#PBS -j eo
#PBS -W group_list=g`id -g`
# ------------------------------
EOF
     ;;
   tropic)
cat >>$RunScript <<EOF
# ------------------------------
#PBS -l ncpus=1
#PBS -l walltime=12:00:00
#PBS -l mem=1gb
#PBS -S /bin/csh
#PBS -V
#PBS -j eo
# ------------------------------
EOF
     ;;
    f01*)
cat >>$RunScript <<EOF
#
# NCCS IBM
#
#@ job_type = parallel
#@ job_name = fvtest
#@ output = \$(job_name).out.\$(jobid)
#@ error =  \$(job_name).err.\$(jobid)
#@ class = poe
#@ tasks_per_node = 1
#@ node = 4
#@ wall_clock_limit = 0:30:00
#@ network.MPI = css0,shared,us
#@ queue
EOF
     ;;
    *)
esac
#
cat >>$RunScript <<EOF
###main() {  
setenv PATH "\${HOME}/bin:\$PATH"     #path to makdep
setenv PATH "/usr/local/bin/:\$PATH"  #path to gmake on phalanx
set LOG=$LOG
set EXELOG=$EXELOG
set Host=$Host
set GPhy=$GPhy
EOF
}

#
#
#
run_fvtest(){
get_headers4RunScript
#
cat >> $RunScript <<EOF

set count=0
set JOBIDS=""  #no space
set JOBIDS2=""
foreach case ($CASES)
   set count=\`expr \$count + 1 \`
   set RunDir=${Home}/../\${case}   
   cd \$RunDir
   gmake cleanall
   if ( "\$count" == 1 ) then
      sed -e "s/NYMDE    = ${Sdatep1}/NYMDE    = ${Edate}/" -e "s/qsub/#qsub/" ./fvgcm.j > ./fv_\${case}.j
   else
      sed -e "s/NYMDE    = ${Sdatep1}/NYMDE    = ${Edate}/"  -e "s/NDAY     = 1/NDAY     = 2/" ./fvgcm.j > ./fv_\${case}.j
   endif
  
   rm -rf d_rst p_rst lsm_rst lsm.rpointer  clm.rpointer clm_rst
   cp /share/fvccm/FVGCM/TEST_DATA/${CASE}/d_rst .
   cp /share/fvccm/FVGCM/TEST_DATA/${CASE}/p_rst .
   if ( "\$GPhy" == 2 ) then
      cp -f /share/fvccm/FVGCM/TEST_DATA/${CASE}/lsm_rst .
      cp -f /share/fvccm/FVGCM/TEST_DATA/${CASE}/lsm.rpointer .
   else
      echo "\$RunDir/clm_rst" >                clm.rpointer 
      cp -f /share/fvccm/FVGCM/TEST_DATA/${CASE}/clm_rst .
   endif

   set mjobid=\`qsub gmake.j\`
   echo "\$mjobid" >> $IDS
   if ( "\$count" == 1 ) then
#    set jobid0=\`qsub -W depend=afterok:\${mjobid} fv_\${case}.j\`
#     set jobid=\`qsub -h -W depend=afterok:\${jobid0} fv_\${case}.j\` #restart
    set jobid0=\`qsub -W depend=afterany:\${mjobid} fv_\${case}.j\`
     echo "\$jobid0" >> $IDS
     set jobid=\`qsub -h -W depend=afterany:\${jobid0} fv_\${case}.j\` #restart
     echo "\$jobid"  >> $IDS
   else
#     set jobid=\`qsub -h -W depend=afterok:\${mjobid} fv_\${case}.j\`
     set jobid=\`qsub -h -W depend=afterany:\${mjobid} fv_\${case}.j\`
     echo "\$jobid"  >> $IDS
   endif
   set JOBIDS="\${JOBIDS}:\${jobid}"
   set JOBIDS2="\${JOBIDS2} \${jobid}"
end
### need to pass  JOBIDS
### chmod +x $RunScript
### qsub $RunScript
### cd $Home
### need to pass  JOBIDS
### go_cmp
### }  #end of run
#
##qsub -W depend=afterany:\${JOBIDS} ${Home}/cmp.csh
cd $Home
set jobid=\`qsub -W depend=afterany\${JOBIDS} ./cmp.csh\`
echo "\$jobid"  >> $IDS
qrls \${JOBIDS2}
EOF
#
go_cmp   #generate cmp.csh before qsub RunSript
#
chmod +x $RunScript
case $Host in
   dycore)
     echo "Running $RunScript "
     #jobid=`qsub $RunScript`  #not in the background
     qsub $RunScript
     ;;
    *)
     echo "Submitting $RunScript "
     jobid=`qsub $RunScript`
     echo "$jobid" >> $IDS
     ;;
esac

}

go_cmp(){
echo "#!/bin/csh -f" > ${Home}/cmp.csh
case $Host in
   dycore)
     ;;
   kalnay|daley|mintz)
cat >>${Home}/cmp.csh <<EOF2
# ------------------------------
#PBS -l ncpus=16
#PBS -l walltime=12:00:00
#PBS -l mem=4gb
#PBS -S /bin/csh
#PBS -V
#PBS -j eo
#PBS -W group_list=g`id -g`
# ------------------------------
EOF2
     ;;
   tropic)
cat >>${Home}/cmp.csh <<EOF2
# ------------------------------
#PBS -l ncpus=1
#PBS -l walltime=1:00:00
#PBS -l mem=4gb
#PBS -S /bin/csh
#PBS -V
#PBS -j eo
#PBS -W group_list=fvccm
# ------------------------------
EOF2
     ;;
    *) 
     ;;
esac
#
cat >>${Home}/cmp.csh <<EOF2
set Phy = $GPhy
set STATUS = \` echo PASSED FAILED\<-- \`
set VerfTag = $TESTDATA
set Host = $Host
switch (\$Host)
   case dycore:
     set OPT="-Wl,-T"
     alias qstat "ps -ef"
     breaksw
   default:
     set OPT=""
     breaksw
endsw
if ( "\$Phy" == 1 ) then
  set str="clm"
else
  set str="lsm"
endif
set VerfDir   = /output/fvgcm/FVGCM/\${VerfTag}/\${str}
set CMPLOG=${Home}/.cmp.log
set CMPLOGTMP=${Home}/.cmp.tmp
rm -f \$CMPLOG \$CMPLOGTMP
cp -f /dev/null \$CMPLOG
set log = \`qstat | grep \$LOGNAME | grep fv_ \`
if ( "\$log" != "" ) then
  echo "Warning: You still have \$log running"
  echo "Date: \`date\`"
endif
set count=0
echo "===================================================="   >>\$CMPLOG
#echo "  Reports of fvtest.j     Tag: ${Ver}     "  >>\$CMPLOG
echo "                 Reports of fvtest.j                "   >>\$CMPLOG
echo "===================================================="   >>\$CMPLOG
foreach case ($CASES)
  set count=\`expr \$count + 1\`
  if ( \$count == 1 )  then
     set S_WORKDIR=$Home/../\$case
     set S_CASE=\$case
     continue
  else
     set WORKDIR=$Home/../\$case 
  endif
# 
# Test 1
#
#  set string = \`echo \$str1   | awk '{printf "%s20",\$1}' \` #not good
  set str1 = "diff \${S_CASE}/d_rst  \${case}/d_rst :"
  set ddate = \`\${WORKDIR}/rst_date \$OPT \${WORKDIR}/d_rst | cut -d" " -f1\`
  if ( "\$ddate" == "${Edate}" ) then
     diff \${S_WORKDIR}/d_rst \${WORKDIR}/d_rst
     set Status = \$status
  else
     set Status = 1
  endif
  set index = \` expr \$Status + 1 \`
  echo "\$str1 \$STATUS[\$index]"  >>\$CMPLOGTMP
#
# Test 2
#
  set str1 = "diff \${S_CASE}/p_rst \${case}/p_rst :"
  set ddate = \`\${WORKDIR}/rst_date \$OPT \${WORKDIR}/d_rst | cut -d" " -f1\`
  if ( "\$ddate" == "${Edate}" ) then
     diff \${S_WORKDIR}/p_rst \${WORKDIR}/p_rst
     set Status = \$status
  else
     set Status = 1
  endif
  set index = \` expr \$Status + 1 \`
  echo "\$str1 \$STATUS[\$index]"  >>\$CMPLOGTMP
end   #end of foreach
#
# Test 3
#
if (-x \${VerfDir}/d_rst  && -x \${VerfDir}/p_rst) then
foreach case ($CASES)
  set WORKDIR=$Home/../\$case 
  set str1 = "diff \${VerfTag}/d_rst \${case}/d_rst :"
  diff \${VerfDir}/d_rst \${WORKDIR}/d_rst
  set Status = \$status
  set index = \` expr \$Status + 1 \`
  echo "\$str1 \$STATUS[\$index]"  >>\$CMPLOGTMP
#
  set str1 = "diff \${VerfTag}/p_rst \${case}/p_rst :"
  diff \${VerfDir}/p_rst \${WORKDIR}/p_rst
  set Status = \$status
  set index = \` expr \$Status + 1 \`
  echo "\$str1 \$STATUS[\$index]"  >>\$CMPLOGTMP
end
endif
#
#
#
awk 'BEGIN{FS=":"}  {printf "%-45s  %s\n",  \$1, \$2}' \$CMPLOGTMP >> \$CMPLOG
echo "                   "  >>\$CMPLOG
echo "Working Dir: $Home"  >>\$CMPLOG
echo "Machine:    \`hostname\` " >> \$CMPLOG
echo "Date:       \`date\`     " >> \$CMPLOG
#/usr/sbin/mailx -s "Reports of fvtest.sh" \$LOGNAME@dao.gsfc.nasa.gov < \$CMPLOG
#/ford1/local/bin/elm -s "Reports of fvtest.sh" \$LOGNAME@dao.gsfc.nasa.gov < \$CMPLOG
elm -s "Reports of fvtest.sh" \$LOGNAME@dao.gsfc.nasa.gov < \$CMPLOG
cat \$CMPLOG
EOF2
chmod +x ${Home}/cmp.csh
}
#
#
#
create_qsub(){
cat >${Home}/qsub <<EOF
#!/bin/csh -f
set idx=\$#argv
set exec=\$argv[\$idx]
echo \$exec
chmod +x \$exec
set filename=\`basename \$exec \`
if ( "\$filename" == "fvtest.j" ) then
  \$exec &
else if ( "\$exec" == "gmake.j" ) then
  \$exec >& /dev/null
else
  \$exec > \${exec}.log
endif
exit 0
EOF
chmod +x ${Home}/qsub

}
#
#
#
export GCase GPhy
GCase=""
GPhy="2"
for arg
do
  case $arg in
    -lsm)
      GPhy=2 ;;
    -clm)
      GPhy=1 ;;
    -case=*)
      GCase=`echo $arg|sed 's/-*case=//'` 
      echo "not implemented"
      exit
      ;;
    -h|-help)
      echo "Usage: fvtest.sh [-lsm] [-clm] [-kill] [-case=b55] "
      exit ;;
    -kill)
      for id in `cat ./.ids`
      do
         echo "qdel $id"
         qdel $id
      done
      exit ;;
     *)
      GPhy=2 ;;
  esac
done
main
