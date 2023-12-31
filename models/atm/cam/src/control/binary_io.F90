#include <misc.h>
#include <params.h>

module binary_io
!----------------------------------------------------------------------- 
! 
! Purpose: wrapper routines for integer, real, SPMD, and non-SPMD binary IO
! 
! Author: 
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8, r4 => shr_kind_r4
   use pmgrid
#if ( defined SPMD )
   use spmd_dyn, only: npes, nlat_p, compute_gsfactors
   use mpishorthand
#endif

CONTAINS

   subroutine wrtout_r8 (iu, arr, numperlat)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Wrapper routine to write restart binary file 
! 
! Author: 
! 
!-----------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in) :: iu                 ! Logical unit
      integer, intent(in) :: numperlat          ! Length of arr

      real(r8) :: arr(numperlat*numlats)        ! Array to be gathered
!
! Local workspace
!
      integer ioerr              ! errorcode

#if ( defined SPMD )
      real(r8), allocatable :: bufres(:) 

      integer :: numrecv(0:npes-1)! number of items to be received
      integer :: displs(0:npes-1) ! displacement array
      integer :: numsend          ! number of bytes to send
      integer :: isiz             ! size of gather array to allocate
      integer :: p                ! process index
      integer :: ier              ! allocation status

      call compute_gsfactors (numperlat, numsend, numrecv, displs)

      if (masterproc) then
         isiz = 0
         do p=0,npes-1
            isiz = isiz + numperlat*nlat_p(p)
         end do
      else
         isiz = 1
      end if
      allocate (bufres(isiz), stat=ier)
      if (ier/=0) call endrun
      
      call mpigatherv (arr, numsend, mpir8, bufres, numrecv, &
                       displs, mpir8, 0, mpicom)

      if (masterproc) then
         write (iu,iostat=ioerr) bufres
         if (ioerr /= 0 ) then
            write (6,*) 'WRTOUT_R8 ioerror ',ioerr,' on i/o unit = ',iu
            call endrun
         end if
         deallocate (bufres)
      endif

#else

      write (iu,iostat=ioerr) arr
      if (ioerr /= 0 ) then
         write (6,*) 'WRTOUT_R8 ioerror ',ioerr,' on i/o unit = ',iu
         call endrun
      end if

#endif

      return
   end subroutine wrtout_r8

!#######################################################################

   subroutine wrtout_r4 (iu, arr, numperlat)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Wrapper routine to write restart binary file 
! 
! Author: 
! 
!-----------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in) :: iu                 ! Logical unit
      integer, intent(in) :: numperlat          ! Length of arr

      real(r4) :: arr(numperlat*numlats)        ! Array to be gathered
!
! Local workspace
!
      integer ioerr              ! errorcode

#if ( defined SPMD )
      real(r4), allocatable :: bufres(:) 

      integer :: numrecv(0:npes-1)! number of items to be received
      integer :: displs(0:npes-1) ! displacement array
      integer :: numsend
      integer :: isiz
      integer :: p

      call compute_gsfactors (numperlat, numsend, numrecv, displs)

      if (masterproc) then
         isiz = 0
         do p=0,npes-1
            isiz = isiz + numperlat*nlat_p(p)
         end do
         allocate (bufres(isiz))
      end if
      
      call mpigatherv (arr, numsend, mpir4, bufres, numrecv, &
                       displs, mpir4, 0, mpicom)

      if (masterproc) then
         write (iu,iostat=ioerr) bufres
         if (ioerr /= 0 ) then
            write (6,*) 'WRTOUT_R4 ioerror ',ioerr,' on i/o unit = ',iu
            call endrun
         end if
         deallocate (bufres)
      endif

#else

      write (iu,iostat=ioerr) arr
      if (ioerr /= 0 ) then
         write (6,*) 'WRTOUT_R4 ioerror ',ioerr,' on i/o unit = ',iu
         call endrun
      end if

#endif

      return
   end subroutine wrtout_r4

!#######################################################################

   subroutine wrtout_int (iu, arr, numperlat)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Wrapper routine to write restart binary file 
! 
! Author: 
! 
!----------------------------------------------------------------------- 
!
! Arguments
!
      integer, intent(in) :: iu                 ! Logical unit
      integer, intent(in) :: numperlat          ! number of values per latitude band
      integer arr(numperlat*numlats)            ! Array to be gathered
!
! Local workspace
!
      integer ioerr              ! errorcode

#if ( defined SPMD )
      integer, allocatable :: bufres(:) 
      integer :: isiz
      integer :: p
      integer :: numsend          ! number of items to be sent
      integer :: displs(0:npes-1) ! displacement array
      integer :: numrecv(0:npes-1)! number of items to be received

      if (masterproc) then
         isiz = 0
         do p=0,npes-1
            isiz = isiz + numperlat*nlat_p(p)
         end do
         allocate (bufres(isiz))
      end if

      call compute_gsfactors (numperlat, numsend, numrecv, displs)
      call mpigatherv (arr, numsend, mpiint, bufres, numrecv, &
                       displs, mpiint, 0, mpicom)

      if (masterproc) then
         write (iu,iostat=ioerr) bufres
         if (ioerr /= 0 ) then
            write (6,*) 'WRTOUT_INT ioerror ',ioerr,' on i/o unit = ',iu
            call endrun
         end if
         deallocate (bufres)
      endif

#else

      write (iu,iostat=ioerr) arr
      if (ioerr /= 0 ) then
         write (6,*) 'WRTOUT_INT ioerror ',ioerr,' on i/o unit = ',iu
         call endrun
      end if
#endif
      return
   end subroutine wrtout_int

!#######################################################################

   subroutine readin_r8 (iu, arr, numperlat)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Wrapper routine to read binary file 
! 
! Author: 
! 
!-----------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in) :: iu                 ! Logical unit
      integer, intent(in) :: numperlat          ! Length of arr

      real(r8) :: arr(numperlat*numlats)        ! Array to be gathered
!
! Local workspace
!
      integer ioerr              ! errorcode

#if ( defined SPMD )
      real(r8), allocatable :: bufres(:) 

      integer :: numrecv          ! number of items to be received
      integer :: displs(0:npes-1) ! displacement array
      integer :: numsend(0:npes-1)
      integer :: isiz
      integer :: p

      call compute_gsfactors (numperlat, numrecv, numsend, displs)

      if (masterproc) then
         isiz = 0
         do p=0,npes-1
            isiz = isiz + numperlat*nlat_p(p)
         end do
         allocate (bufres(isiz))
         read (iu,iostat=ioerr) bufres
         if (ioerr /= 0 ) then
            write (6,*) 'READIN_R8 ioerror ',ioerr,' on i/o unit = ',iu
            call endrun
         end if
      end if
      
      call mpiscatterv (bufres, numsend, displs, mpir8, arr, &
                        numrecv, mpir8, 0, mpicom)

      if (masterproc) then
         deallocate (bufres)
      endif

#else

      read (iu,iostat=ioerr) arr
      if (ioerr /= 0 ) then
         write (6,*) 'READIN_R8 ioerror ',ioerr,' on i/o unit = ',iu
         call endrun
      end if

#endif

      return
   end subroutine readin_r8

!#######################################################################

   subroutine readin_r4 (iu, arr, numperlat)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Wrapper routine to read binary file 
! 
! Author: 
! 
!-----------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in) :: iu                 ! Logical unit
      integer, intent(in) :: numperlat          ! Length of arr

      real(r4) :: arr(numperlat*numlats)        ! Array to be gathered
!
! Local workspace
!
      integer ioerr              ! errorcode

#if ( defined SPMD )
      real(r4), allocatable :: bufres(:) 

      integer :: numrecv          ! number of items to be received
      integer :: displs(0:npes-1) ! displacement array
      integer :: numsend(0:npes-1)
      integer :: isiz
      integer :: p

      call compute_gsfactors (numperlat, numrecv, numsend, displs)

      if (masterproc) then
         isiz = 0
         do p=0,npes-1
            isiz = isiz + numperlat*nlat_p(p)
         end do
         allocate (bufres(isiz))
         read (iu,iostat=ioerr) bufres
         if (ioerr /= 0 ) then
            write (6,*) 'READIN_R4 ioerror ',ioerr,' on i/o unit = ',iu
            call endrun
         end if
      end if
      
      call mpiscatterv (bufres, numsend, displs, mpir4, arr, &
                        numrecv, mpir4, 0, mpicom)

      if (masterproc) then
         deallocate (bufres)
      endif

#else

      read (iu,iostat=ioerr) arr
      if (ioerr /= 0 ) then
         write (6,*) 'READIN_R4 ioerror ',ioerr,' on i/o unit = ',iu
         call endrun
      end if

#endif

      return
   end subroutine readin_r4

!#######################################################################

   subroutine readin_int (iu, arr, numperlat)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Wrapper routine to write restart binary file 
! 
! Author: 
! 
!-----------------------------------------------------------------------
!
! Arguments
!
      integer iu                 ! Logical unit
      integer numperlat          ! number of values per latitude band
      integer arr(numperlat*numlats)        ! Array to be gathered
!
! Local workspace
!
      integer ioerr              ! errorcode

#if ( defined SPMD )
      integer, allocatable :: bufres(:) 
      integer isiz
      integer :: p
      integer :: numrecv          ! number of items to be received
      integer :: displs(0:npes-1) ! displacement array
      integer :: numsend(0:npes-1)! number of items to be sent

      call compute_gsfactors (numperlat, numrecv, numsend, displs)

      if (masterproc) then
         isiz = 0
         do p=0,npes-1
            isiz = isiz + numperlat*nlat_p(p)
         end do
         allocate (bufres(isiz))
         read (iu,iostat=ioerr) bufres
         if (ioerr /= 0 ) then
            write (6,*) 'READIN_INT ioerror ',ioerr,' on i/o unit = ',iu
            call endrun
         end if
      end if

      call mpiscatterv (bufres, numsend, displs, mpiint, arr, &
                        numrecv, mpiint, 0, mpicom)

      if (masterproc) then
         deallocate (bufres)
      endif

#else

      read (iu,iostat=ioerr) arr
      if (ioerr /= 0 ) then
         write (6,*) 'READIN_INT ioerror ',ioerr,' on i/o unit = ',iu
         call endrun
      end if
#endif
      return
   end subroutine readin_int

end module binary_io
