#!/bin/csh
#depends on results of build_them happening successfully.
#setenv REGNUM 128
setenv REGNUM 32
set thisdir = `pwd`
#foreach theregamt ( 0.01 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 )
#foreach theregamt ( 0.01 0.1 0.15 0.2 0.25 0.3 0.35 0.4 )
#foreach theregamt ( 0.005 0.0025 )
foreach theregamt ( 0.01 0.25 )
  set xx = "${REGNUM}x${theregamt}"
  set exename = "nncamreg_${xx}"
  set rundir=$SCRATCH/regularization_suite_03/$xx
  mkdir -p $rundir
  cp ${thisdir}/*.i.* $rundir
  cp -r ${thisdir}/keras_matrices $rundir
  cat ${thisdir}/atm_in.template | sed "s@REGNUM@${REGNUM}@g" | sed "s@REGAMT@${theregamt}@g" > $rundir/atm_in
  cp ${thisdir}/exe/$exename $rundir
  cat ${thisdir}/run.slurm.template | sed "s@RUNDIR@${rundir}@g" | sed "s@REGNUM@${REGNUM}@g" | sed "s@REGAMT@${theregamt}@g" | sed "s@EXENAME@${exename}@g" > ${rundir}/run.slurm.$xx
  cd ${rundir}
#  sbatch ./run.slurm.$xx
end
