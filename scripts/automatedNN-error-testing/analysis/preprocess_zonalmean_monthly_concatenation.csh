#!/bin/csh
module load nco
set outputdir="/scratch/08110/tg874091/automatedNN-error-testing"
set controloutputdir="/scratch/08110/tg874091/smoketestNN/run"
set resultsdir=$SCRATCH/test # HEY jerry update this path to wherever you like.
echo on
# HEY Jerry When you are done
# scp all zonalmean_*.nc in the resultdir to wherever you are doing the matlab analysis.
mkdir -p $resultsdir
foreach varname ( T Q PS )
  foreach k ( `seq 1 108` )
#  foreach k ( 3 3 )
    set modlabel = `printf "%05d" $k`
    echo "$modlabel"
  # Step 1 -- extract single variable and concatenate across monthly raw output files into temporary file...
    set rcatfile = "${resultsdir}/${modlabel}_${varname}_rcat.nc"
    set h0filelist = "$outputdir/run_${modlabel}/*.cam2.h0.*.nc"
    echo "$h0filelist"
    if ( ! -e $rcatfile ) then
      ncrcat -D 2 -v hyam,hybm,hyai,hybi,$varname $h0filelist -o ${rcatfile}
      # note the hy* stuff is needed downstream for calculation of the pressure thickness field from the surface pressure.
    endif
       # Step 2 --- further average over the longitude dimension to create the zonal mean monthly climatology (which is what we will analyze)
    set zmfile = "${resultsdir}/zonalmean_${modlabel}_${varname}_rcat.nc"
    if ( ! -e $zmfile ) then
      if ( -e $rcatfile ) then
        ncwa -a lon $rcatfile -o $zmfile
      endif
    endif
  end #k, loop over NN tests
  # same steps for the control run:
  echo " ------- PROCESSING CONTROL RUN -------"
  set rcatfile = "$resultsdir/control_${varname}_rcat.nc"
  if ( ! -e $rcatfile ) then
    ncrcat -v hyam,hybm,hyai,hybi,$varname $controloutputdir/*.cam2.h0.*.nc -o ${rcatfile}
  endif
  set zmfile = "${resultsdir}/zonalmean_control_${varname}_rcat.nc"
  if ( ! -e $zmfile ) then
    ncwa -a lon $rcatfile -o $zmfile
  endif
end

