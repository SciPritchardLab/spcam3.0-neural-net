#! /usr/bin/csh -f

##-----------------------------------------------------------------------
## Hopper
##------------
##

## netCDF stuff
setenv INC_NETCDF   /opt/cray/netcdf/3.6.2/netcdf-pgi/include
setenv LIB_NETCDF   /opt/cray/netcdf/3.6.2/netcdf-pgi/lib

## $wrkdir is a working directory where the model will be built and run.
## $blddir is the directory where model will be compiled.
## $rundir is the directory where the model will be run.
## $cfgdir is the directory containing the CAM configuration scripts.
set wrkdir       = `pwd`/obj
set blddir       = $wrkdir
set rundir       = `pwd`/run
set cfgdir       = `pwd`/bld

## Ensure that run and build directories exist
mkdir -p $rundir                || echo "cannot create $rundir" && exit 1
mkdir -p $blddir                || echo "cannot create $blddir" && exit 1

## If an executable doesn't exist, build one.
if ( ! -x $blddir/cam ) then
    cd $blddir                  || echo "cd $blddir failed" && exit 1
    $cfgdir/configure_mmf -fc ftn -spmd -nosmp -dyn sld -res 64x128 -pcols 8 -nlev 30 -cam_exedir $rundir -mpi_inc $MPICH_DIR/include -mpi_lib $MPICH_DIR/lib || echo "configure failed" && exit 1
    echo "building CAM in $blddir ..."
    rm -f Depends
    gmake >&! MAKE.out      || echo "CAM build failed: see $blddir/MAKE.out" && exit 1
endif

exit 0
