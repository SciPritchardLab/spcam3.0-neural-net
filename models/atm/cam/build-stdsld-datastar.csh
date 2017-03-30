#! /usr/bin/csh -f

#-----------------------------------------------------------------------
## IBM
##------------
##

## netCDF stuff
setenv INC_NETCDF   /usr/local/apps64/netcdf/netcdf-3.6.1/include
setenv LIB_NETCDF   /usr/local/apps64/netcdf/netcdf-3.6.1/lib

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
    $cfgdir/configure -dyn sld -res 64x128 -pcols 8 -nlev 26 -cam_exedir $rundir || echo "configure failed" && exit 1
    echo "building CAM in $blddir ..."
    rm -f Depends
    gmake -j16 #>&! MAKE.out      || echo "CAM build failed: see $blddir/MAKE.out" && exit 1
endif

exit 0
