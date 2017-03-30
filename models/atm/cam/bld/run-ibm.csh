#! /usr/bin/csh -f

#-----------------------------------------------------------------------
## IBM
##------------
##
## This is an example script to build and run the default CAM configuration
## on an IBM SP.  The default configuration is T42L26, Eulerian dynamics,
## CLM2 land model, and CSIM4 ice model.
##
## Setting LoadLeveler options for batch queue submission.
## @ class is the queue in which to run the job.
##     To see a list of available queues, type "llclass" interactively.
## @ node is the number of nodes. @tasks_per_node should be set to 1 because
##     of the hybrid OpenMP/MPI configuration of CAM.  The number of nodes should
##     be a power of 2, up to a max of 16 for T42.
## @ output and @error are the names of file written to the directory from
##     which the script is submitted containing STDOUT and STDERR respectively
## @ job_type = parallel declares that multiple nodes will be used.
## @ network.MPI: Has to do with network connection between nodes.  Best to leave alone.
## @ node_usage = not_shared acquires dedicated access to nodes for the job.
## @ queue tells load leveler to submit the job

#@ class          = com_reg
#@ node           = 2
#@ tasks_per_node = 1
#@ output         = out.$(jobid)
#@ error          = out.$(jobid)
#@ job_type       = parallel
#@ network.MPI    = csss,not_shared,us
#@ node_usage     = not_shared
#@ queue

## POE Environment.  Set these for interactive jobs.  They're ignored by LoadLeveler
## MP_NODES is the number of nodes.  The number chosen should be a power of 2, up to a max of 16 for T42.
setenv MP_NODES 2
setenv MP_TASKS_PER_NODE 1
setenv MP_EUILIB us
setenv MP_RMPOOL 1
# TH:  bug fix suggested by Brian Eaton  1/24/03
unsetenv MP_PROCS

## suggestion from Jim Edwards to reintroduce XLSMPOPTS on 11/13/03
setenv XLSMPOPTS "stack=256000000"
setenv AIXTHREAD_SCOPE S
setenv MALLOCMULTIHEAP true
setenv OMP_DYNAMIC false
## Do our best to get sufficient stack memory
limit stacksize unlimited

## netCDF stuff
setenv INC_NETCDF /usr/local/include
setenv LIB_NETCDF /usr/local/lib64/r4i4

## ROOT OF CAM DISTRIBUTION - probably needs to be customized.
## Contains the source code for the CAM distribution.
## (the root directory contains the subdirectory "models")
set camroot      = /fs/cgd/data0/$LOGNAME/cam2_0

## ROOT OF CAM DATA DISTRIBUTION - needs to be customized unless running at NCAR.
## Contains the initial and boundary data for the CAM distribution.
## (the root directory contains the subdirectories "atm" and "lnd")
setenv CSMDATA    /fs/cgd/csm/inputdata

## Default namelist settings:
## $case is the case identifier for this run. It will be placed in the namelist.
## $runtype is the run type: initial, restart, or branch.
## $nelapse is the number of timesteps to integrate, or number of days if negative.
set case         = camrun
set runtype      = initial
set nelapse      = -1

## $wrkdir is a working directory where the model will be built and run.
## $blddir is the directory where model will be compiled.
## $rundir is the directory where the model will be run.
## $cfgdir is the directory containing the CAM configuration scripts.
set wrkdir       = /ptmp/$LOGNAME
set blddir       = $wrkdir/$case/bld
set rundir       = $wrkdir/$case
set cfgdir       = $camroot/models/atm/cam/bld

## Ensure that run and build directories exist
mkdir -p $rundir                || echo "cannot create $rundir" && exit 1
mkdir -p $blddir                || echo "cannot create $blddir" && exit 1

## If an executable doesn't exist, build one.
if ( ! -x $blddir/cam ) then
    cd $blddir                  || echo "cd $blddir failed" && exit 1
    $cfgdir/configure           || echo "configure failed" && exit 1
    echo "building CAM in $blddir ..."
    rm -f Depends
    gmake -j4 >&! MAKE.out      || echo "CAM build failed: see $blddir/MAKE.out" && exit 1
endif

## Create the namelist
cd $blddir                      || echo "cd $blddir failed" && exit 1
$cfgdir/build-namelist -s -case $case -runtype $runtype -o $rundir/namelist \
 -namelist "&camexp nelapse=$nelapse /"  || echo "build-namelist failed" && exit 1

## Run CAM
cd $rundir                      || echo "cd $rundir failed" && exit 1
echo "running CAM in $rundir"
poe $blddir/cam < namelist      || echo "CAM run failed" && exit 1

exit 0
