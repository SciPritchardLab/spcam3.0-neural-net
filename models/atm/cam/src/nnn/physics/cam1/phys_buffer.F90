#include <misc.h>
#include <params.h>

module phys_buffer

!----------------------------------------------------------------------- 
! 
! Purpose: 
!   Implement a physics buffer to hold arrays that must persist
!   across timesteps or between calls to different physics packages.
!
!   Current implementation only supports one buffer.
!
! Author: B. Eaton
! 
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use infnan,       only: nan
   use pmgrid,       only: masterproc
   use ppgrid,       only: pcols, begchunk, endchunk
   use phys_grid,    only: physgrid_set, write_field_from_chunk, read_chunk_from_field
   use prognostics,  only: ptimelevels
   use string_utils, only: to_upper
#ifdef SPMD
   use mpishorthand, only: mpicom, mpiint
#endif

   implicit none
   private
   save

! Public methods 

   public ::&
      pbuf_init,          &! initialize physics buffer
      pbuf_add,           &! add field to physics buffer
      pbuf_get_fld_idx,   &! get index of specified field in the physics buffer
      pbuf_old_tim_idx,   &! return the index for the oldest time
      pbuf_next_tim_idx,  &! return the index for the next time
      pbuf_update_tim_idx,&! update the index for the oldest time
      pbuf_allocate,      &! allocate memory for physics buffer fields
      pbuf_deallocate,    &! deallocate memory for physics buffer fields
      pbuf_write_restart, &! write physics buffer to restart file
      pbuf_read_restart    ! read physics buffer from restart file

! Public types and data

   type, public :: pbuf_fld
      character(len=16)                       :: name
      character(len=16)                       :: scope
      integer                                 :: fdim, mdim, ldim
      real(r8), pointer, dimension(:,:,:,:,:) :: fld_ptr
   end type pbuf_fld

   integer, public, parameter :: pbuf_size_max=100
   type(pbuf_fld), public, dimension(pbuf_size_max) :: pbuf

   integer, public :: pbuf_times  ! number of time levels in physics buffer (dycore dependent)

! Private module data

   integer :: pbuf_size = 0
   integer :: old_time_idx = 1

!=========================================================================================
contains
!=========================================================================================

subroutine pbuf_init()

! Initialize physics buffer.

   implicit none
!-----------------------------------------------------------------------------------------

   pbuf_times = ptimelevels - 1

end subroutine pbuf_init

!=========================================================================================

subroutine pbuf_add(name, scope, fdim, mdim, ldim, index)

! Add a field to the physics buffer

   implicit none

   character(len=*), intent(in)  :: name   ! field name 
   character(len=*), intent(in)  :: scope  ! field scope, either 'global' or 'physpkg'
   integer,          intent(in)  :: fdim   ! first generic field dimension
   integer,          intent(in)  :: mdim   ! middle generic field dimension
   integer,          intent(in)  :: ldim   ! last generic field dimension
   integer,          intent(out) :: index  ! index in the physics buffer

! Local variables
   character(len=*), parameter :: sub = 'pbuf_add'
   integer :: i
!-----------------------------------------------------------------------------------------

   if ( pbuf_size >= pbuf_size_max ) then
      write(6,*)sub,': max number of physics buffer fields exceeded.'
      write(6,*)'  Increase pbuf_size_max parameter in phys_buffer module.'
      call endrun
   end if

   do i = 1, pbuf_size
      if ( pbuf(i)%name == name ) then
         write(6,*)sub,': ERROR: field name ', name, ' is already in use.'
         call endrun
      end if
   end do

   pbuf_size = pbuf_size + 1
   index = pbuf_size
   pbuf(index)%name = name
   pbuf(index)%scope = scope
   pbuf(index)%fdim = fdim
   pbuf(index)%mdim = mdim
   pbuf(index)%ldim = ldim

end subroutine pbuf_add

!=========================================================================================

function pbuf_get_fld_idx(name)

! Get index of specified field in the physics buffer.  String matching is case insensitive.
! Call endrun if name not found

   implicit none

   character(len=*), intent(in)  :: name   ! field name 

! Return value
   integer :: pbuf_get_fld_idx

! Local variables
   integer :: i
!-----------------------------------------------------------------------------------------

   do i = 1, pbuf_size
      if ( to_upper(pbuf(i)%name) == to_upper(name) ) then
         pbuf_get_fld_idx = i
         return
      end if
   end do

   write(6,*)'ERROR in pbuf_get_fld_idx: index not found for ',name
   call endrun

end function pbuf_get_fld_idx

!=========================================================================================

function pbuf_old_tim_idx()

! Return index of oldest time sample in the physics buffer.

   implicit none

! Return value
   integer :: pbuf_old_tim_idx
!-----------------------------------------------------------------------------------------

   pbuf_old_tim_idx = old_time_idx

end function pbuf_old_tim_idx

!=========================================================================================

function pbuf_next_tim_idx(idx)

! Return index of next time sample in the physics buffer.

   implicit none

   integer, intent(in) :: idx

! Return value
   integer :: pbuf_next_tim_idx
!-----------------------------------------------------------------------------------------

   pbuf_next_tim_idx = mod(idx, pbuf_times) + 1

end function pbuf_next_tim_idx

!=========================================================================================

subroutine pbuf_update_tim_idx()

! Update index of old time sample in the physics buffer.

   implicit none
!-----------------------------------------------------------------------------------------

   old_time_idx = mod(old_time_idx, pbuf_times) + 1

end subroutine pbuf_update_tim_idx

!=========================================================================================

subroutine pbuf_allocate(scope)

! Allocate storage for fields in the physics buffer with the specified scope.
!
! N.B. This routine must be called after phys_grid_init because that's
!      where begchunk and endchunk are set

   implicit none

   character(len=*), intent(in)  :: scope

! Local variables
   character(len=*), parameter :: sub = 'pbuf_allocate'
   integer :: i, fdim, mdim, ldim, istat
!-----------------------------------------------------------------------------------------

   if ( .not. physgrid_set ) then
      write(6,*)sub,': ERROR: called before physics grid initialized.'
      call endrun
   end if

   do i = 1, pbuf_size
      if ( pbuf(i)%scope == scope ) then
         fdim = pbuf(i)%fdim
         mdim = pbuf(i)%mdim
         ldim = pbuf(i)%ldim
         allocate(pbuf(i)%fld_ptr(fdim,pcols,mdim,begchunk:endchunk,ldim), stat=istat)
         if ( istat /= 0 ) then
            write(6,*)sub,': ERROR: allocate failed for ',pbuf(i)%name
            call endrun
         end if
         pbuf(i)%fld_ptr = nan
      end if
   end do

end subroutine pbuf_allocate

!=========================================================================================

subroutine pbuf_deallocate(scope)

! Deallocate storage for fields in the physics buffer with the specified scope.

   implicit none

   character(len=*), intent(in)  :: scope

! Local variables
   character(len=*), parameter :: sub = 'pbuf_deallocate'
   integer :: i, fdim, mdim, ldim
!-----------------------------------------------------------------------------------------

   do i = 1, pbuf_size
      if ( pbuf(i)%scope == scope ) then
         if (associated(pbuf(i)%fld_ptr)) then
            deallocate(pbuf(i)%fld_ptr)
         else
            write(6,*)sub,': ERROR: ',pbuf(i)%name, ' is not allocated'
            call endrun
         end if
      end if
   end do

end subroutine pbuf_deallocate

!=========================================================================================

subroutine pbuf_write_restart(iu)

! write physics buffer to restart file

   implicit none

   integer, intent(in) :: iu  ! Fortran unit number

! Local variables
   character(len=*), parameter :: sub = 'pbuf_write_restart'
   integer :: i, ioerr
!-----------------------------------------------------------------------------------------

   if (masterproc) then
      write (iu,iostat=ioerr) old_time_idx
      if (ioerr /= 0 ) then
         write (6,*) sub,' ioerror ', ioerr,' on i/o unit = ',iu
         call endrun
      end if
   endif

   do i = 1, pbuf_size
      if ( pbuf(i)%scope == 'global' ) then
         call write_field_from_chunk(iu, pbuf(i)%fdim, pbuf(i)%mdim, pbuf(i)%ldim, pbuf(i)%fld_ptr)
      end if
   end do

end subroutine pbuf_write_restart

!=========================================================================================

subroutine pbuf_read_restart(iu)

! write physics buffer to restart file

   implicit none

   integer, intent(in) :: iu  ! Fortran unit number

! Local variables
   character(len=*), parameter :: sub = 'pbuf_read_restart'
   integer :: i, ioerr
!-----------------------------------------------------------------------------------------

   if (masterproc) then
      read (iu,iostat=ioerr) old_time_idx
      if (ioerr /= 0 ) then
         write (6,*) sub,' ioerror ', ioerr,' on i/o unit = ',iu
         call endrun
      end if
   endif

#if ( defined SPMD ) 
   call mpibcast(old_time_idx, 1, mpiint, 0, mpicom)
#endif

   do i = 1, pbuf_size
      if ( pbuf(i)%scope == 'global' ) then
         call read_chunk_from_field(iu, pbuf(i)%fdim, pbuf(i)%mdim, pbuf(i)%ldim, pbuf(i)%fld_ptr)
      end if
   end do

end subroutine pbuf_read_restart

!=========================================================================================

end module phys_buffer
