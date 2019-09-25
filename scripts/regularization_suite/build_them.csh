#!/bin/csh
#setenv REGNUM 128
setenv REGNUM 32
set thisdir = `pwd`
set srcdir = "/home1/00993/tg802402/repositories/spcam3.0-neural-net/models/atm/cam/src/physics/cam1"
# this scripts assumes that we have already built the model once here:
set objdir = "/home1/00993/tg802402/nncambld/obj"
set tmpexedir = "/home1/00993/tg802402/nncambld/run"
# ... such that we can just tweak one source file and repeatedly do a partil re-build:

#foreach theregamt ( 0.01 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 )
#foreach theregamt ( 0.005 0.0025 )
foreach theregamt ( 0.01 0.25 )
  set xx = "${REGNUM}x${theregamt}"
  set exefile = "$thisdir/exe/nncamreg_${xx}"
  if ( ! -e $exefile ) then
    cat ${thisdir}/cloudbrain.F90.template | sed "s@XXX@${REGNUM}@g" | sed "s@YYY@${theregamt}@g" > $srcdir/cloudbrain.F90
    cd $objdir 
    gmake
    cp $tmpexedir/cam $exefile
  endif
end
