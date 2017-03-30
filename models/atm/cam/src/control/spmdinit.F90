#include <misc.h>
#include <params.h>

subroutine spmdinit
!----------------------------------------------------------------------- 
! 
! Purpose: MPI initialization routine:  
! 
! Method: get number of cpus, processes, tids, etc
!         dynamics and physics decompositions are set up later
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------

#if ( defined SPMD )
   use mpishorthand, only: mpiint, mpichar, mpilog, mpipk, mpir8, mpir4, &
                           mpicom, mpi_max_processor_name, mpi_integer, &
                           mpi_character, mpi_logical, mpi_real8, mpi_real4, &
                           mpi_packed, mpi_comm_world
#if ( defined IRIX64 )
   use mpishorthand, only: mpi_thread_multiple
#endif
   use spmd_dyn, only:  npes, nsmps, proc_smp_map, spmdinit_dyn
#endif
   use pmgrid, only: plat, masterproc, iam
   use dycore, only: dycore_is

   implicit none

#if ( defined SPMD )
!
! Local workspace
!
   integer i,j               ! indices
   integer npthreads         ! thread status
   integer ier               ! return error status    
   integer, allocatable :: length(:)  ! length of name
   logical is_sgi, done
   
   character*(mpi_max_processor_name), allocatable :: proc_name(:) ! returned processor name
   character*(mpi_max_processor_name), allocatable :: smp_name(:)  ! SMP name

   logical mpi_running       ! returned value indicates if MPI_INIT has been called
!---------------------------------------------------------------------------
!
! Initialize mpi
!
#if (! defined COUP_CSM)
   call mpi_initialized (mpi_running, ier)
!
! If SGI and fv_dycore, initialize MPI with mpi_init_thread in order for
! mod_comm MPI-2 threaded calls to work (Mirin, Sawyer)
!
#if ( defined IRIX64 )
   is_sgi = .true.
#else
   is_sgi = .false.
#endif
   if (.not.mpi_running) then
      if (dycore_is ('LR') .and. is_sgi) then
#if ( defined IRIX64 )
         print *, 'Initializing MPI for FV dycore on SGI'
         call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE, npthreads, ier)
         call MPI_QUERY_THREAD(npthreads, ier)
         if (npthreads /= MPI_THREAD_MULTIPLE) then
            write(*,*) 'did not provide MPI_THREAD_MULTIPLE. ', &
                  'Change to MPI_THREAD_MULTIPLE with MPI_INIT_THREAD ', &
                  'for multi-threading MPI2'
         endif
#endif
      else
         call mpi_init (ier)       
      endif
   endif
#endif

!
! Set mpishorthand variables.  Need to set as variables rather than parameters since
! some MPI implementations set values for MPI tags at run time
!
   mpiint  = mpi_integer
   mpichar = mpi_character
   mpilog  = mpi_logical
   mpir4   = mpi_real4
   mpir8   = mpi_real8
   mpipk   = mpi_packed
#if (defined COUP_CSM)
!  mpicom now set in cam.F90 for COUP_CSM
!  using cpl_interface_init()
#else
   call mpi_comm_dup(mpi_comm_world, mpicom, ier)
#endif
!
! Get my id  
!
   call mpi_comm_rank (mpicom, iam, ier) 
   if (iam == 0) then 
      masterproc = .true.
   else
      masterproc = .false.
   end if
!
! Get number of processors
!
   call mpi_comm_size (mpicom, npes, ier)

   allocate ( length(0:npes-1) )
   allocate ( proc_name(0:npes-1) )
   proc_name(:) = ' '

!
! Get processor names and send to root. "1" is the msg tag
!
   call mpi_get_processor_name (proc_name(iam), length(iam), ier)

   if (masterproc) then
      do i=1,npes-1
         call mpirecv (proc_name(i), mpi_max_processor_name, mpichar, i, 1, mpicom)
      end do
      write(6,*) npes, 'pes participating in computation'
      write(6,*) '-----------------------------------'
      write(6,*) 'NODE#  NAME'
      do i=0,npes-1
         write(6,'(i3,2x,a)') i,trim(proc_name(i))
      end do
   else
      call mpisend (proc_name(iam), mpi_max_processor_name, mpichar, 0, 1, mpicom)
   end if
!
! Identify SMP nodes and process/SMP mapping.
! (Assume that processor names are SMP node names on SMP clusters.)
!
   allocate ( proc_smp_map(0:npes-1) )
   if (masterproc) then
      allocate ( smp_name(0:npes-1) )
      smp_name(:) = ' '
      proc_smp_map(:) = -1
!
      nsmps = 1
      smp_name(0) = trim(proc_name(0))
      proc_smp_map(0) = 0
!
      do i=1,npes-1
         j = 0
         done = .false.
         do while ((.not. done) .and. (j < nsmps))
            if (smp_name(j) .eq. trim(proc_name(i))) then
               proc_smp_map(i) = j
               done = .true.
            endif
            j = j + 1
         enddo
         if (.not. done) then
            smp_name(nsmps) = trim(proc_name(i))
            proc_smp_map(i) = nsmps
            nsmps = nsmps + 1
         endif
      enddo
      deallocate(smp_name)
   endif
   call mpibcast(nsmps, 1, mpiint, 0, mpicom)
   call mpibcast(proc_smp_map, npes, mpiint, 0, mpicom)
!
   deallocate(length)
   deallocate(proc_name)
#endif   
   return
end subroutine spmdinit

