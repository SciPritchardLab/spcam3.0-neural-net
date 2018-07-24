#!/bin/csh
set rawdir = "/scratch/00993/tg802402/cloudbrain_ctrl_aquaplanet"
foreach varname ( PRECT FLUT ) #SPDT SPDQ QRL QRS )
  set scratchfile = "./scratch/${varname}_rcat.nc"
  ncrcat -D 2 -v $varname $rawdir/*.h1.000?-??-21-00000.nc -o $scratchfile  &
end
wait

