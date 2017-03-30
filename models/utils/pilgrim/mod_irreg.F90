!------------------------------------------------------------------------------
!BOP
!
! !MODULE: mod_irreg --- Additions to mod_comm for irregular communication
      module mod_irreg
!
! !USES:
#include "pilgrim.h"
#if !defined( STAND_ALONE )
#include "misc.h"
#endif
#if defined ( SPMD )
      use mod_comm, only : mp_r8, mp_r4, mp_i4, window,                        &
                     r8_win, i4_win, ga_r8_s, ga_r8_r, ga_i4_s, ga_i4_r,       &
                     gid, numpro, numcpu, pkgs_per_pro, commglobal, Status,    &
                     Win_Open, Win_Close, Win_Finalize,                        &
                     Ga_RecvInit_r8, Ga_RecvInit_i4

#if !defined(USE_MLP) && !defined(MPI2)
      use mod_comm, only : sqest, rqest, nsend, nread
#endif

#if defined( TIMING )
      use timingModule, only : timing_on, timing_off
#endif
#if defined( STAND_ALONE )
#define r8 selected_real_kind(12)
#define r4 selected_real_kind( 6)
#define i8 selected_int_kind(13)
#define i4 selected_int_kind( 6)
#else
      use shr_kind_mod, only : r8 => shr_kind_r8, r4 => shr_kind_r4,  &
                               i8 => shr_kind_i8, i4 => shr_kind_i4
#endif

      implicit none

#if !defined ( USE_MLP )
#include "mpif.h"
#endif
 
! !PUBLIC MEMBER FUNCTIONS:
      public mp_sendirr, mp_recvirr, mp_sendirr_i4, mp_recvirr_i4
      public blockdescriptor

!------------------------------------------------------------------------------
!  type declaration for describing an arbitrary number of contiguous chunks
!------------------------------------------------------------------------------
      type blockdescriptor
#if defined( USE_MPI_TYPES )
         integer              :: type               ! Ptr to MPI derived type
#else
         integer, pointer     :: displacements(:)   ! Offsets in local segment
         integer, pointer     :: blocksizes(:)      ! Block sizes to transfer
         integer              :: partneroffset      ! Aggregated partner offset
#endif
      end type blockdescriptor

! !DESCRIPTION:
!
! !REVISION HISTORY:
!    2003.06.04   Sawyer  Created from CAM mod_comm.F90
!    2003.06.05   Sawyer  Mirin bug fix -- exploit full buffer, not just half
!    2003.06.24   Sawyer  Integrated USE_MPI_TYPES from parutilitiesmodule
!
!EOP
!------------------------------------------------------------------------------
      INTEGER, SAVE :: InHandle(MAX_PAX, MAX_TRF)
      INTEGER, SAVE :: OutHandle(MAX_PAX, MAX_TRF)
      INTEGER, SAVE :: BegTrf = 0  ! Ongoing overlapped begintransfer #
      INTEGER, SAVE :: EndTrf = 0  ! Ongoing overlapped endtransfer #

      contains

!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_sendirr --- Write r8 contiguous chunks to global array
!
! !INTERFACE:
      subroutine mp_sendirr ( q, send_bl, recv_bl, qout )
 
! !INPUT PARAMETERS:
      type(blockdescriptor), intent(in)  :: send_bl(:) ! send blocks
      type(blockdescriptor), intent(in)  :: recv_bl(:) ! receive blocks
      real(r8), intent(in) :: q(*)                     ! local data segment

! !OUTPUT PARAMETERS:
      real(r8), intent(out) :: qout(*)                 ! local output segment
!
! !DESCRIPTION:
!     Copy a number of contiguous chunks into GA.  This fundamental
!     routine forms a basis for higher level primitives
!
! !REVISION HISTORY: 
!    02.08.13   Sawyer      Creation
!    02.11.06   Mirin       Optimizations
!    03.03.03   Sawyer      Use partneroffset
!    03.06.24   Sawyer      Integrated USE_MPI_TYPES; added qout
!
! !BUGS:
!    MLP not yet running
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer ipe, qsize, offset, blocksize, nchunks, offset_s, offset_r, ierr

!
!     initialize window
#if defined ( USE_MLP )
#include "mlp_ptr.h"
#endif

#if defined( TIMING )
      call timing_on('COMM_TOTAL')
#endif

#if defined( USE_MPI_TYPES )
!
! Increment the ongoing transfer number
      BegTrf = MOD(BegTrf,MAX_TRF) + 1
!
! MPI: Irecv over all processes
!
      do ipe=1, numpro

!
! Receive the buffers with MPI_Irecv. Non-blocking
!
        OutHandle(ipe,BegTrf) = MPI_REQUEST_NULL
        if ( recv_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
          call mpi_irecv( qout, 1, recv_bl(ipe)%type, ipe-1, ipe-1,     &
                          commglobal, OutHandle(ipe,BegTrf), ierr )
        endif
      enddo

!
! MPI: Isend over all processes
!
      do ipe=1, numpro

!
! Send the individual buffers with non-blocking sends
!
        InHandle(ipe,BegTrf) = MPI_REQUEST_NULL
        if ( send_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
          call mpi_isend( q, 1, send_bl(ipe)%type, ipe-1, gid,        &
                          commglobal, InHandle(ipe,BegTrf), ierr )
        endif
      enddo
#else

      call Win_Open(r8_win)
      blocksize = r8_win%size / numpro   ! size alloted to this PE
      offset = gid*blocksize
      offset_s = 0
      offset_r = 0

      do ipe=1, numpro
#if !defined(USE_MLP) && !defined(MPI2)
         r8_win%size_r = SUM(recv_bl(ipe)%blocksizes)
         if (r8_win%size_r .ne. 0) then
            r8_win%offset_r = offset_r
            offset_r = offset_r + r8_win%size_r
            if (offset_r .gt. r8_win%size) then
              print *, "Fatal mp_sendirr: receive window out of space - exiting"
              call exit(1)
            endif
            r8_win%src = ipe-1
            call Ga_RecvInit_r8(r8_win, ga_r8_r)
         endif
#endif
         qsize = SUM( send_bl(ipe)%blocksizes )
         if (qsize .ne. 0) then
            r8_win%dest = ipe-1
#if defined(USE_MLP)
            r8_win%offset_s = offset
#else
            r8_win%offset_s = offset_s
            offset_s = offset_s + qsize
            if (offset_s .gt. r8_win%size) then
              print *, "Fatal mp_sendirr: send window out of space - exiting"
              call exit(1)
            endif
#if defined(MPI2)
            r8_win%offset_r = send_bl(ipe)%partneroffset
#endif
#endif
            nchunks = size( send_bl(ipe)%displacements )
            call Ga_PutContig_r8( q, r8_win, nchunks, &
              send_bl(ipe)%displacements, send_bl(ipe)%blocksizes, ga_r8_s )
         endif
         offset = offset + qsize
      enddo
#endif

#if defined( TIMING )
      call timing_off('COMM_TOTAL')
#endif

      end subroutine mp_sendirr
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_recvirr --- Read r8 contiguous chunks to global array
!
! !INTERFACE:
      subroutine mp_recvirr ( qout, recv_bl )
 
! !INPUT PARAMETERS:
      type(blockdescriptor), intent(in)  :: recv_bl(:) ! Global Array Window
! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent(inout) :: qout(*)               ! local data segment
!
! !DESCRIPTION:
!
!     Complete transfer of a generalized region.
!
! !REVISION HISTORY:
!    02.08.15   Sawyer      Creation
!    02.11.06   Mirin       Optimizations
!    03.03.03   Sawyer      Now using packed arrays for MPI2
!
!EOP
!------------------------------------------------------------------------------
!BOC
      integer :: ipe, blocksize, nchunks, offset_r
#if defined( USE_MPI_TYPES)
      integer Ierr
      integer InStats(numpro*MPI_STATUS_SIZE)
      integer OutStats(numpro*MPI_STATUS_SIZE)
#endif


#if defined ( USE_MLP ) && defined ( LAHEY )
#include "mlp_ptr.h"
#endif

#if defined( TIMING )
      call timing_on('COMM_TOTAL')
#endif

#if defined(USE_MPI_TYPES)
      EndTrf = MOD(EndTrf,MAX_TRF) + 1
      CALL MPI_WAITALL( numpro, InHandle(:,EndTrf), InStats, Ierr )
      CALL MPI_WAITALL( numpro, OutHandle(:,EndTrf), OutStats, Ierr )
#else
      call Win_Close(r8_win)

      blocksize = r8_win%size / numpro   ! size alloted to this PE
      offset_r = 0
      do ipe=1, numpro
         r8_win%size_r = SUM(recv_bl(ipe)%blocksizes)
         if (r8_win%size_r .ne. 0) then
#if defined(USE_MLP)
            r8_win%offset_r = offset_r + XXXXXXXXXXXXXXXXXX !  NOT YET TESTED!!
#else
            r8_win%offset_r = offset_r
            offset_r = offset_r + r8_win%size_r
            if (offset_r .gt. r8_win%size) then
              print *, "Fatal mp_recvirr: receive window out of space - exiting"
              call exit(1)
            endif
#endif
            nchunks = SIZE( recv_bl(ipe)%displacements )
            call Ga_GetContig_r8 ( qout, r8_win, nchunks,                      &
            recv_bl(ipe)%displacements, recv_bl(ipe)%blocksizes, ga_r8_r )
         endif
      enddo

      call Win_Finalize(r8_win)
#endif

#if defined( TIMING )
      call timing_off('COMM_TOTAL')
#endif

!EOC
      end subroutine mp_recvirr
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_sendirr_i4 --- Write i4 contiguous chunks to global array
!
! !INTERFACE:
      subroutine mp_sendirr_i4 ( q, send_bl, recv_bl, qout )
 
! !INPUT PARAMETERS:
      type(blockdescriptor), intent(in)  :: send_bl(:) ! send blocks
      type(blockdescriptor), intent(in)  :: recv_bl(:) ! receive blocks
      integer(i4), intent(in) :: q(*)                  ! local data segment

! !OUTPUT PARAMETERS:
      integer(i4), intent(out) :: qout(*)              ! local output segment
!
! !DESCRIPTION:
!     Copy a number of contiguous chunks into GA.  This fundamental
!     routine forms a basis for higher level primitives
!
! !REVISION HISTORY: 
!    02.08.13   Sawyer      Creation
!    02.11.06   Mirin       Optimizations
!    03.03.03   Sawyer      Use partneroffset
!    03.06.24   Sawyer      Integrated USE_MPI_TYPES; added qout
!
! !BUGS:
!    MLP not yet running
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer ipe, qsize, offset, blocksize, nchunks, offset_s, offset_r, ierr

!
!     initialize window
#if defined ( USE_MLP )
#include "mlp_ptr.h"
#endif

#if defined( TIMING )
      call timing_on('COMM_TOTAL')
#endif

#if defined(USE_MPI_TYPES)
!
! Increment the ongoing transfer number
      BegTrf = MOD(BegTrf,MAX_TRF) + 1
!
! MPI: Irecv over all processes
!
      do ipe=1, numpro

!
! Receive the buffers with MPI_Irecv. Non-blocking
!
        OutHandle(ipe,BegTrf) = MPI_REQUEST_NULL
        if ( recv_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
          call mpi_irecv( qout, 1, recv_bl(ipe)%type, ipe-1, ipe-1,     &
                          commglobal, OutHandle(ipe,BegTrf), ierr )
        endif
      enddo

!
! MPI: Isend over all processes
!
      do ipe=1, numpro

!
! Send the individual buffers with non-blocking sends
!
        InHandle(ipe,BegTrf) = MPI_REQUEST_NULL
        if ( send_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
          call mpi_isend( q, 1, send_bl(ipe)%type, ipe-1, gid,        &
                          commglobal, InHandle(ipe,BegTrf), ierr )
        endif
      enddo
#else
      call Win_Open(i4_win)
      blocksize = i4_win%size / numpro   ! size alloted to this PE
      offset = gid*blocksize
      offset_s = 0
      offset_r = 0

      do ipe=1, numpro
#if !defined(USE_MLP) && !defined(MPI2)
         i4_win%size_r = SUM(recv_bl(ipe)%blocksizes)
         if (i4_win%size_r .ne. 0) then
            i4_win%offset_r = offset_r
            offset_r = offset_r + i4_win%size_r
            if (offset_r .gt. i4_win%size) then
              print *, "Fatal mp_sendirr: receive window out of space - exiting"
              call exit(1)
            endif
            i4_win%src = ipe-1
            call Ga_RecvInit_i4(i4_win, ga_i4_r)
         endif
#endif
         qsize = SUM( send_bl(ipe)%blocksizes )
         if (qsize .ne. 0) then
            i4_win%dest = ipe-1
#if defined(USE_MLP)
            i4_win%offset_s = offset
#else
            i4_win%offset_s = offset_s
            offset_s = offset_s + qsize
            if (offset_s .gt. i4_win%size) then
              print *, "Fatal mp_sendirr: send window out of space - exiting"
              call exit(1)
            endif
#if defined(MPI2)
            i4_win%offset_r = send_bl(ipe)%partneroffset
#endif
#endif
            nchunks = size( send_bl(ipe)%displacements )
            call Ga_PutContig_i4( q, i4_win, nchunks, &
              send_bl(ipe)%displacements, send_bl(ipe)%blocksizes, ga_i4_s )
         endif
         offset = offset + qsize
      enddo
#endif

#if defined( TIMING )
      call timing_off('COMM_TOTAL')
#endif

      end subroutine mp_sendirr_i4
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_recvirr_i4 --- Read i4 contiguous chunks to global array
!
! !INTERFACE:
      subroutine mp_recvirr_i4 ( qout, recv_bl )
 
! !INPUT PARAMETERS:
      type(blockdescriptor), intent(in)  :: recv_bl(:) ! Global Array Window
! !INPUT/OUTPUT PARAMETERS:
      integer(i4), intent(inout) :: qout(*)               ! local data segment
!
! !DESCRIPTION:
!
!     Complete transfer of a generalized region.
!
! !REVISION HISTORY:
!    02.08.15   Sawyer      Creation
!    02.11.06   Mirin       Optimizations
!    03.03.03   Sawyer      Now using packed arrays for MPI2
!
!EOP
!------------------------------------------------------------------------------
!BOC
      integer :: ipe, blocksize, nchunks, offset_r
#if defined( USE_MPI_TYPES)
      integer Ierr
      integer InStats(numpro*MPI_STATUS_SIZE)
      integer OutStats(numpro*MPI_STATUS_SIZE)
#endif

#if defined ( USE_MLP ) && defined ( LAHEY )
#include "mlp_ptr.h"
#endif

#if defined( TIMING )
      call timing_on('COMM_TOTAL')
#endif

#if defined(USE_MPI_TYPES)
      EndTrf = MOD(EndTrf,MAX_TRF) + 1
      CALL MPI_WAITALL( numpro, InHandle(:,EndTrf), InStats, Ierr )
      CALL MPI_WAITALL( numpro, OutHandle(:,EndTrf), OutStats, Ierr )
#else
      call Win_Close(i4_win)

      blocksize = i4_win%size / numpro     ! size alloted to this PE
      offset_r = 0
      do ipe=1, numpro
         i4_win%size_r = SUM(recv_bl(ipe)%blocksizes)
         if (i4_win%size_r .ne. 0) then
#if defined(USE_MLP)
            i4_win%offset_r = offset_r + XXXXXXXXXXXXXXXXXX !  NOT YET TESTED!!
#else
            i4_win%offset_r = offset_r
            offset_r = offset_r + i4_win%size_r
            if (offset_r .gt. i4_win%size) then
              print *, "Fatal mp_recvirr: receive window out of space - exiting"
              call exit(1)
            endif
#endif
            nchunks = SIZE( recv_bl(ipe)%displacements )
            call Ga_GetContig_i4 ( qout, i4_win, nchunks,                      &
            recv_bl(ipe)%displacements, recv_bl(ipe)%blocksizes, ga_i4_r )
         endif
      enddo

      call Win_Finalize(i4_win)
#endif

#if defined( TIMING )
      call timing_off('COMM_TOTAL')
#endif

!EOC
      end subroutine mp_recvirr_i4
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_PutContig_r8 --- Write r8 contiguous chunks to global array
!
! !INTERFACE:
      subroutine Ga_PutContig_r8 ( q, win, nchunks, displ, sizes, ga )
 
! !INPUT PARAMETERS:
      type(window), intent(in)  :: win           ! Global Array Window
      integer, intent(in)  :: nchunks            ! number chunks
      integer, intent(in)  :: displ(nchunks)     ! starting displacements
      integer, intent(in)  :: sizes(nchunks)     ! chunk lengths
      real(r8), intent(in) :: q(*)               ! local data segment
! !INPUT/OUTPUT PARAMETERS:
#if defined(USE_MLP)
      real(r8), intent(inout):: ga(win%size*max_call)
#else
      real(r8), intent(inout):: ga(win%size)
#endif
!
! !DESCRIPTION:
!     Copy a number of contiguous chunks into GA.  This fundamental
!     routine forms a basis for higher level primitives
!
! !REVISION HISTORY: 
!    02.08.13   Sawyer        Creation
!    03.03.03   Sawyer        Added OpenMP
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, j, qsize, mydisp, send_tag, ierror, p, mysize, tmpsize

      mydisp = win%offset_s
      qsize = 0
      do j = 1, nchunks
         do i = displ(j)+1, displ(j)+sizes(j)
            qsize = qsize+1
            ga(mydisp+qsize) = q(i)
         enddo
      enddo

#if !defined(USE_MLP)
#if defined(MPI2)
      if ( win%offset_r+qsize > win%size ) then
         print *, "Fatal PutContig: window out of space - exiting"
         print *, 'gid win%offset_r qsize win%size = ', gid,      &
                   win%offset_r, qsize, win%size
         call exit(1)
      endif
      tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
!$omp parallel do private(p,mysize,ierror)
      do p=1,MIN(pkgs_per_pro,qsize)
        mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
        call MPI_PUT(ga(win%offset_s+(p-1)*tmpsize+1), mysize, mp_r8,     &
                     win%dest, win%offset_r+(p-1)*tmpsize, mysize, mp_r8, &
                     win%id, ierror)
      enddo
#else
      send_tag = gid
      nsend = nsend + 1
      call MPI_ISEND(ga(mydisp+1), qsize, mp_r8, win%dest, &
                     send_tag, commglobal, sqest(nsend), ierror)
#endif
#endif
      end subroutine Ga_PutContig_r8
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_GetContig_r8 --- Read from real*8 4d global array
!
! !INTERFACE:
      subroutine Ga_GetContig_r8 ( q, win, nchunks, displ, sizes, ga )
! !INPUT PARAMETERS:
      type(window), intent(in)  :: win           ! Global Array Window
      integer, intent(in)  :: nchunks            ! number chunks
      integer, intent(in)  :: displ(nchunks)     ! starting displacements
      integer, intent(in)  :: sizes(nchunks)     ! chunk lengths
#if defined(USE_MLP)
      real(r8), intent(in) :: ga(win%size*max_call)
#else
      real(r8), intent(in) :: ga(win%size)
#endif
! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent(inout) :: q(*)
!
! !DESCRIPTION:
!
!     Read contiguous chunks from a real*8 global array.  This fundamental
!     routine forms a basis for higher level primitives
!
! !REVISION HISTORY: 
!    02.08.13   Sawyer        Creation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, j, mydisp, qsize, ierror

#if !defined(USE_MLP) && !defined(MPI2)
      nread = nread + 1
      call MPI_WAIT(rqest(nread), Status, ierror)
#endif

      mydisp = win%offset_r
      qsize = 0
      do j = 1, nchunks
         do i = displ(j)+1, displ(j)+sizes(j)
            qsize = qsize+1
            q(i) = ga(mydisp+qsize)
         enddo
      enddo
!!!      print *, "GetContig:gid", gid, "mydisp", mydisp, "size", qsize, "nchunks", nchunks
#if defined(USE_MLP)
      if ( qsize > win%size/(2.*numpro) ) then
         print *, "Fatal GetContig: window out of space - exiting"
         call exit(1)
      endif
#endif
!EOC
      end subroutine Ga_GetContig_r8
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_PutContig_i4 --- Write i4 contiguous chunks to global array
!
! !INTERFACE:
      subroutine Ga_PutContig_i4 ( q, win, nchunks, displ, sizes, ga )
 
! !INPUT PARAMETERS:
      type(window), intent(in)  :: win           ! Global Array Window
      integer, intent(in)  :: nchunks            ! number chunks
      integer, intent(in)  :: displ(nchunks)     ! starting displacements
      integer, intent(in)  :: sizes(nchunks)     ! chunk lengths
      integer(i4), intent(in) :: q(*)               ! local data segment
! !INPUT/OUTPUT PARAMETERS:
#if defined(USE_MLP)
      integer(i4), intent(inout):: ga(win%size*max_call)
#else
      integer(i4), intent(inout):: ga(win%size)
#endif
!
! !DESCRIPTION:
!     Copy a number of contiguous chunks into GA.  This fundamental
!     routine forms a basis for higher level primitives
!
! !REVISION HISTORY: 
!    02.08.13   Sawyer        Creation
!    03.03.03   Sawyer        Added OpenMP
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, j, qsize, mydisp, send_tag, ierror, p, mysize, tmpsize

      mydisp = win%offset_s
      qsize = 0
      do j = 1, nchunks
         do i = displ(j)+1, displ(j)+sizes(j)
            qsize = qsize+1
            ga(mydisp+qsize) = q(i)
         enddo
      enddo

#if !defined(USE_MLP)
#if defined(MPI2)
      if ( win%offset_r+qsize > win%size ) then
         print *, "Fatal PutContig: window out of space - exiting"
         call exit(1)
      endif
      tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
!$omp parallel do private(p,mysize,ierror)
      do p=1,MIN(pkgs_per_pro,qsize)
        mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
        call MPI_PUT(ga(win%offset_s+(p-1)*tmpsize+1), mysize, mp_i4,     &
                     win%dest, win%offset_r+(p-1)*tmpsize, mysize, mp_i4, &
                     win%id, ierror)
      enddo
#else
      send_tag = gid
      nsend = nsend + 1
      call MPI_ISEND(ga(mydisp+1), qsize, mp_i4, win%dest, &
                     send_tag, commglobal, sqest(nsend), ierror)
#endif
#endif
      end subroutine Ga_PutContig_i4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_GetContig_i4 --- Read from i4 4d global array
!
! !INTERFACE:
      subroutine Ga_GetContig_i4 ( q, win, nchunks, displ, sizes, ga )
! !INPUT PARAMETERS:
      type(window), intent(in)  :: win           ! Global Array Window
      integer, intent(in)  :: nchunks            ! number chunks
      integer, intent(in)  :: displ(nchunks)     ! starting displacements
      integer, intent(in)  :: sizes(nchunks)     ! chunk lengths
#if defined(USE_MLP)
      integer(i4), intent(in) :: ga(win%size*max_call)
#else
      integer(i4), intent(in) :: ga(win%size)
#endif
! !INPUT/OUTPUT PARAMETERS:
      integer(i4), intent(inout) :: q(*)
!
! !DESCRIPTION:
!
!     Read contiguous chunks from a int*4 global array.  This fundamental
!     routine forms a basis for higher level primitives
!
! !REVISION HISTORY: 
!    02.08.13   Sawyer        Creation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, j, mydisp, qsize, ierror

#if !defined(USE_MLP) && !defined(MPI2)
      nread = nread + 1
      call MPI_WAIT(rqest(nread), Status, ierror)
#endif

      mydisp = win%offset_r
      qsize = 0
      do j = 1, nchunks
         do i = displ(j)+1, displ(j)+sizes(j)
            qsize = qsize+1
            q(i) = ga(mydisp+qsize)
         enddo
      enddo
!!!      print *, "GetContig:gid", gid, "mydisp", mydisp, "size", qsize, "nchunks", nchunks
#if defined(USE_MLP)
      if ( qsize > win%size/(2.*numpro) ) then
         print *, "Fatal GetContig: window out of space - exiting"
         call exit(1)
      endif
#endif
!EOC
      end subroutine Ga_GetContig_i4
!------------------------------------------------------------------------------

#endif
      end module mod_irreg

