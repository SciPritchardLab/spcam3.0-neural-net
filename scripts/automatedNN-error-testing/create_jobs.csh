#!/bin/csh
set pathOttmodels = "/work/08110/tg874091/stampede2/all_models_fkb_ott_et_al"
#foreach k ( `seq 1 108` )
foreach k ( `seq 1 10` )
  set modlabel = `printf "%05d" $k`
  set rundir = "run_${modlabel}"
  mkdir -p $rundir
  cp -r baseline/* $rundir 
# Copy in Jordan's NN to overwrite keras_matrices/model.txt
  cp $pathOttmodels/$modlabel.txt $rundir/keras_matrices/model.txt
  cat atm_in.template | sed "s@XXX@$modlabel@g" > $rundir/atm_in
  cat run.template | sed "s@XXX@$modlabel@g" > $rundir/run.csh
  sbatch $rundir/run.csh
end
