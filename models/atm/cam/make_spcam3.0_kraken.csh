#! /usr/bin/csh -f
# Mike Pritchard, University of Washington
# mspritch@uw.edu
# Created Feb. 10, 2012
# Build and run the model on NICS Kraken, a Cray XT5

# Note default compilers on login are already PGI but accesssed through a 
# wrapper called "ftn" (passed with -fc flag to configure below, and activated through a new clause in the models/atm/cam/bld/Makefile.in file

module load netcdf/3.6.2
setenv INC_NETCDF   $NETCDF_DIR/include
setenv LIB_NETCDF   $NETCDF_DIR/lib
setenv CSMDATA /lustre/scratch/mikeprit/csmdata_cam3.0

#set camroot = $HOME/src/spcam3.0/cam3_sp
set camroot = $HOME/src/spcam3.0_mikemods/cam3_sp
set cfgdir  = $camroot/models/atm/cam/bld
echo $cfgdir

set wrkdir       = `pwd`/obj
set blddir       = $wrkdir
set rundir       = `pwd`/run

## Ensure that run and build directories exist
mkdir -p $rundir                || echo "cannot create $rundir" && exit 1
mkdir -p $blddir                || echo "cannot create $blddir" && exit 1

# If an executable doesn't exist, build one.
if ( ! -x $blddir/cam ) then
    cd $blddir                  || echo "cd $blddir failed" && exit 1
    $cfgdir/configure_mmf -spmd -smp -fc ftn -dyn sld -res 64x128 -pcols 8 -nlev 30 -cam_exedir $rundir -mpi_inc $MPICH_DIR/include -mpi_lib $MPICH_DIR/lib || echo "configure failed" && exit 1
    echo "building CAM in $blddir ..."
    rm -f Depends
    gmake -j 32 #>&! MAKE.out      || echo "CAM build failed: see $blddir/MAKE.out" && exit 1
endif

cp $blddir/cam $rundir
$cfgdir/build-namelist -o $rundir/atm_in -config $blddir/config_cache.xml -csmdata $CSMDATA -namelist "&camexp ncdata='$CSMDATA/atm/cam/inic/gaus/cami_0000-09-01_64x128_L30_c031210.nc' /"

exit 0

# WARNING the time dimension of history files created is EMPTY unless you
# also do a #define HDEBUG in history.F90 (weird behavior)

# To run the model on Kraken interactively on 64 MPI processes (max for T42 spectral resolution) do this:

#qsub -I -l size=32,walltime=01:00:00 -A TG-XXXXXXX
# aprun -n 32 ./cam < atm_in > logfile &
# tail -f logfile

