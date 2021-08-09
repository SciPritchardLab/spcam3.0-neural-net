#! /usr/bin/csh -f
module load netcdf
# Use Intel MPI libraries
set rpath="$HOME/repos/spcam3.0-neural-net"
setenv camroot $rpath/models/atm/cam
setenv esmfroot $rpath/models/utils/esmf/build/linux_intel
echo $camroot
#exit(0)
cp $camroot/bld/Makefile.stampede $camroot/bld/Makefile
# Note I had to install my own version of netcdf3.6.3 to be old enough to play nice with spcam3.
setenv INC_NETCDF   /home1/00993/tg802402/include 
setenv LIB_NETCDF   /home1/00993/tg802402/lib

setenv MPICH_DIR $MPICH_HOME # on stampede, this is env variable for impi home after module set right.

# user override if desired. Expectation is you will use this scripts in a local dir, where obj created
set wrkdir       = $SCRATCH/smoketestNN
set blddir       = $wrkdir/obj
set rundir       = $wrkdir/run
set cfgdir       = $camroot/bld

## Ensure that run and build directories exist
mkdir -p $rundir                || echo "cannot create $rundir" && exit 1
mkdir -p $blddir                || echo "cannot create $blddir" && exit 1

## If an executable doesn't exist, build one.
if ( ! -x $blddir/cam ) then
    cd $blddir                  || echo "cd $blddir failed" && exit 1

# for SP control run activate this version:
#    $cfgdir/configure_mmf -fc mpif90 -cc cc -spmd -smp -dyn sld -res 64x128 -pcols 8 -nlev 30 -cam_exedir $rundir -mpi_inc $MPICH_DIR/intel64/include -mpi_lib $MPICH_DIR/intel64/lib || echo "configure failed" && exit 1
# for NN run activate this version:
    $cfgdir/configure_mmf -fflags "-DSPFLUXBYPASS -DCLOUDBRAIN -DNEURALLIB -DBRAINDEBUG -DHDEBUG" -fc mpif90 -cc cc -spmd -smp -dyn sld -res 64x128 -pcols 8 -nlev 30 -cam_exedir $rundir -mpi_inc $MPICH_DIR/intel64/include -mpi_lib $MPICH_DIR/intel64/lib || echo "configure failed" && exit 1
# Use this non-SP build script when using CLOUDBRAIN to avoid stomping on state_save compiler messages:
    echo "building CAM in $blddir ..."
    rm -f Depends
    gmake -j 8 >&! Make.out #-j 8 #>&! MAKE.out      || echo "CAM build failed: see $blddir/MAKE.out" && exit 1
endif

exit 0
