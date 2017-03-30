#! /usr/bin/csh -f

# Use Intel MPI libraries
module switch mvapich2 impi
setenv camroot ${HOME}/GitHub/spcam3.0-neural-net/models/atm/cam
setenv esmfroot ${HOME}/GitHub/spcam3.0-neural-net/models/utils/esmf/build/linux_intel
cp $camroot/bld/Makefile.stampede $camroot/bld/Makefile
setenv INC_NETCDF   /opt/apps/intel13/netcdf/3.6.3/x86_64/include
setenv LIB_NETCDF   /opt/apps/intel13/netcdf/3.6.3/x86_64/lib

setenv MPICH_DIR $MPICH_HOME # on stampede, this is env variable for impi home after module set right.

# user override if desired. Expectation is you will use this scripts in a local dir, where obj created
set wrkdir       = `pwd`/obj
set blddir       = $wrkdir
set rundir       = `pwd`/run
set cfgdir       = $camroot/bld

## Ensure that run and build directories exist
mkdir -p $rundir                || echo "cannot create $rundir" && exit 1
mkdir -p $blddir                || echo "cannot create $blddir" && exit 1

## If an executable doesn't exist, build one.
if ( ! -x $blddir/cam ) then
    cd $blddir                  || echo "cd $blddir failed" && exit 1
    $cfgdir/configure_mmf -fc mpif90 -cc cc -spmd -smp -dyn sld -res 64x128 -pcols 8 -nlev 30 -cam_exedir $rundir -mpi_inc $MPICH_DIR/intel64/include -mpi_lib $MPICH_DIR/intel64/lib || echo "configure failed" && exit 1
    echo "building CAM in $blddir ..."
    rm -f Depends
    gmake -j 8 #>&! MAKE.out      || echo "CAM build failed: see $blddir/MAKE.out" && exit 1
endif

exit 0
