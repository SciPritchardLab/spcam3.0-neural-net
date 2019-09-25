#!/bin/csh
mkdir -p ./results
set regn = 128
module load nco
foreach regamt ( 0.0025 0.005 0.01 0.05 0.15 0.25 0.35 0.01 0.1 0.2 0.3 0.4  )
 set xx = ${regn}x${regamt}
 set rundir = $SCRATCH/regularization_suite_03/$xx
 set resultsfile = "./results/TPHYSTND_lev25_rcat_${xx}.nc"
 if ( ! -e $resultsfile ) then
   ncrcat -D 2 -d lev,880.,910. $rundir/*.h1.*.nc -o $resultsfile &
 endif
end
wait


