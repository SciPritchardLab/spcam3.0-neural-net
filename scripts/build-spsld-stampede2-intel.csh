#! /usr/bin/csh -f
module load netcdf
# Use Intel MPI libraries
set rpath="$HOME/repositories/spcam3.0-neural-net_causalcoupler"
setenv camroot $rpath/models/atm/cam
setenv esmfroot $rpath/models/utils/esmf/build/linux_intel
echo $camroot
cp $camroot/bld/Makefile.stampede $camroot/bld/Makefile
# Note I had to install my own version of netcdf3.6.3 to be old enough to play nice with spcam3.
setenv INC_NETCDF   $HOME/include 
setenv LIB_NETCDF   $HOME/lib
setenv MPICH_DIR $MPICH_HOME # on stampede, this is env variable for impi home after module set right.

# user override if desired. Expectation is you will use this scripts in a local dir, where obj created
set casenm       = SingleNNsCAM3.0_test1.4
set nnmodels     = dummy_singleNNs_FKB_renamed_v1.1
set wrkdir       = $HOME/CAM3.0-builds/$casenm
set blddir       = $wrkdir/obj
set rundir       = $HOME/CAM3.0-builds/$casenm/run
set cfgdir       = $camroot/bld
set srundir      = /scratch/08098/tg873976/$casenm
set initfile     = spinup_AndKua_aqua_SPCAM3.0.cam2.i.0000-12-02-00000.nc
set namelist     = test_atm_in

## Ensure that run and build directories exist
mkdir -p $rundir                || echo "cannot create $rundir" && exit 1
mkdir -p $blddir                || echo "cannot create $blddir" && exit 1
mkdir -p $srundir               || echo "cannot create $srundir" && exit 1
ln -s $rpath/run_dir/models/$nnmodels $srundir/models
ln -s $rpath/run_dir/$initfile $srundir/$initfile
ln -s $rpath/run_dir/$namelist $srundir/atm_in
ln -s $rundir/cam $srundir/cam

## If an executable doesn't exist, build one.
if ( ! -x $blddir/cam ) then
    cd $blddir                  || echo "cd $blddir failed" && exit 1

    # NNCAM
    $cfgdir/configure_mmf -fflags "-DDEEP -DCLOUDBRAIN" -fc mpif90 -cc cc -spmd -smp -dyn sld -res 64x128 -pcols 8 -nlev 30 -cam_exedir $rundir -mpi_inc $MPICH_DIR/intel64/include -mpi_lib $MPICH_DIR/intel64/lib || echo "configure failed" && exit 1
    
    # SPCAM
#    $cfgdir/configure_mmf -fc mpif90 -cc cc -spmd -smp -dyn sld -res 64x128 -pcols 8 -nlev 30 -cam_exedir $rundir -mpi_inc $MPICH_DIR/intel64/include -mpi_lib $MPICH_DIR/intel64/lib || echo "configure failed" && exit 1

    # CAM
#    $cfgdir/configure -fc mpif90 -cc cc -spmd -smp -dyn sld -res 64x128 -pcols 8 -nlev 30 -cam_exedir $rundir -mpi_inc $MPICH_DIR/intel64/include -mpi_lib $MPICH_DIR/intel64/lib || echo "configure failed" && exit 1

    echo "building CAM in $blddir ..."
    rm -f Depends
    gmake -j 8 >&! Make.out #-j 8 #>&! MAKE.out      || echo "CAM build failed: see $blddir/MAKE.out" && exit 1
endif

exit 0
