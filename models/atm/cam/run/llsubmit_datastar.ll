#!/usr/bin/ksh
#@environment = COPY_ALL;\
#AIXTHREAD_SCOPE=S;\
#MP_ADAPTER_USE=dedicated;\
#MP_CPU_USE=unique;\
#MP_CSS_INTERRUPT=no;\
#MP_EAGER_LIMIT=64K;\
#MP_EUIDEVELOP=min;\
#MP_LABELIO=yes;\
#MP_POLLING_INTERVAL=100000;\
#MP_PULSE=0;\
#MP_SHARED_MEMORY=yes;\
#MP_SINGLE_THREAD=no;\
#RT_GRQ=ON;\
#SPINLOOPTIME=0;\
#YIELDLOOPTIME=0
#@account_no = cst102
#@class = normal
#@node = 64
#@tasks_per_node = 1
#@wall_clock_limit = 18:00:00
#@node_usage = not_shared
#@network.MPI = sn_all, shared, US
#@job_type = parallel
#@job_name= job.$(jobid)
#@output = LL_out.$(jobid)
#@error = LL_err.$(jobid)
#@notification = always
##@notify_user = mark@atmos.colostate.edu
#@initialdir = /gpfs/bransonm/mmf_cam3.0/cam3_sp/models/atm/cam/run 
#@queue
export MEMORY_AFFINITY=MCM
export MP_TASK_AFFINITY=MCM
export MALLOCMULTIHEAP=heaps:8
export OMP_NUM_THREADS=8
poe cam < namelist

