#include <misc.h>
#include <params.h>

module restart_physics

   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid
   use phys_grid,   only: read_chunk_from_field, write_field_from_chunk, &
                          gather_chunk_to_field, get_ncols_p
   use pmgrid,      only: masterproc
   use prognostics, only: ptimelevels, n3
   use buffer
   use radae,       only: abstot_3d, absnxt_3d, emstot_3d, initialize_radbuffer
   use comsrf
   use ioFileMod
#if ( defined COUP_CSM )
   use ccsm_msg, only: initialize_ccsm_msg, write_restart_ccsm, read_restart_ccsm
#endif
#ifdef QRLDAMP
   use qrl_anncycle, only: qrl_dailymean, qrl_dailymean_buffer, qrl_anncycle_int
#endif

   implicit none

   private
!
! Public interfaces
!
   public write_restart_physics    ! Write the physics restart info out
   public read_restart_physics     ! Read the physics restart info in
   public get_abs_restart_filepath ! Get the name of the restart filepath
!
! Private data
!
   character(len=256) :: pname  ! Full abs-ems restart filepath
!
! Filename specifier for restart abs-ems file
! (%c = caseid, $y = year, $m = month, $d = day, $s = seconds in day, %t = tape number)
!
   character(len=256) :: rafilename_spec = '%c.cam2.ra.%y-%m-%d-%s'   ! abs-ems restart

CONTAINS

   subroutine write_restart_physics (nrg, nrg2)
      use filenames, only: mss_irt, mss_wpass, get_archivedir, interpret_filename_spec
      use phys_buffer, only: pbuf_write_restart
! for nlend and aeres
#include <comctl.h>
!
! Input arguments
!
      integer :: nrg
      integer :: nrg2
!
! Local workspace
!
      real(r8) tmpfield(pcols,begchunk:endchunk)
      real(r8) tmpfield3d(pcols,plevmx,begchunk:endchunk)
      integer i                 ! loop index
      integer n3tmp             ! timestep index
      character(len=256) fname  ! abs-ems restart filename
      integer ioerr             ! I/O status
      integer  :: ncol          ! number of vertical columns
!
! Buffer module variables
!
      call write_field_from_chunk(nrg,1,1,1,pblht)
      call write_field_from_chunk(nrg,1,1,1,tpert)
      call write_field_from_chunk(nrg,1,pver,1,qrs)
      call write_field_from_chunk(nrg,1,pver,1,qrl)
      call write_field_from_chunk(nrg,1,pcnst+pnats,1,qpert)

! Physics buffer
      call pbuf_write_restart(nrg)
!
! Comsrf module variables
!
#if (! defined COUP_CSM)
      call write_field_from_chunk(nrg,1,1,1,fsnt)
#endif
      call write_field_from_chunk(nrg,1,1,1,fsns)
#if (! defined COUP_CSM)
      call write_field_from_chunk(nrg,1,1,1,flnt)
      call write_field_from_chunk(nrg,1,1,1,flns)
#endif
      do i=begchunk,endchunk
	tmpfield(:,i) = srfflx_state2d(i)%asdir(:)
      end do
      call write_field_from_chunk(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
	tmpfield(:,i) = srfflx_state2d(i)%asdif(:)
      end do
      call write_field_from_chunk(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
	tmpfield(:,i) = srfflx_state2d(i)%aldir(:)
      end do
      call write_field_from_chunk(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
	tmpfield(:,i) = srfflx_state2d(i)%aldif(:)
      end do
      call write_field_from_chunk(nrg,1,1,1,tmpfield)
#if (! defined COUP_CSM)
      call write_field_from_chunk(nrg,1,1,1,asdirice)
      call write_field_from_chunk(nrg,1,1,1,asdifice)
      call write_field_from_chunk(nrg,1,1,1,aldirice)
      call write_field_from_chunk(nrg,1,1,1,aldifice)
      call write_field_from_chunk(nrg,1,1,1,asdirocn)
      call write_field_from_chunk(nrg,1,1,1,asdifocn)
      call write_field_from_chunk(nrg,1,1,1,aldirocn)
      call write_field_from_chunk(nrg,1,1,1,aldifocn)
      call write_field_from_chunk(nrg,1,1,1,asdirlnd)
      call write_field_from_chunk(nrg,1,1,1,asdiflnd)
      call write_field_from_chunk(nrg,1,1,1,aldirlnd)
      call write_field_from_chunk(nrg,1,1,1,aldiflnd)
      call write_field_from_chunk(nrg,1,1,1,lwuplnd)
      call write_field_from_chunk(nrg,1,1,1,lwupocn)
      call write_field_from_chunk(nrg,1,1,1,lwupice)
      call write_field_from_chunk(nrg,1,1,1,tsice)
#endif

      do i=begchunk,endchunk
	tmpfield(:,i) = srfflx_state2d(i)%lwup(:)
      end do
      call write_field_from_chunk(nrg,1,1,1,tmpfield)
      call write_field_from_chunk(nrg,1,1,1,landfrac)
#if ( defined COUP_SOM )
      call write_field_from_chunk(nrg,1,1,1,aice)
#endif
      call write_field_from_chunk(nrg,1,1,1,ocnfrac)
      call write_field_from_chunk(nrg,1,1,1,icefrac)
      call write_field_from_chunk(nrg,1,1,1,landm)
      call write_field_from_chunk(nrg,1,1,1,sgh)
      do i=begchunk,endchunk
	tmpfield(:,i) = srfflx_state2d(i)%ts(:)
      end do
      call write_field_from_chunk(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
	tmpfield(:,i) = srfflx_state2d(i)%sst(:)
      end do
      call write_field_from_chunk(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
	tmpfield3d(:,:,i) = surface_state2d(i)%tssub(:,:)
      end do
      call write_field_from_chunk(nrg,1,plevmx,1,tmpfield3d)
      call write_field_from_chunk(nrg,1,1,1,sicthk)
      call write_field_from_chunk(nrg,1,1,1,snowhland)
#if (! defined COUP_CSM)
      call write_field_from_chunk(nrg,1,1,1,snowhice)
#endif
      do i=begchunk,endchunk
       ncol = get_ncols_p(i)
	tmpfield(:ncol,i) = surface_state2d(i)%flwds(:ncol)
      end do
      call write_field_from_chunk(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
       ncol = get_ncols_p(i)
	tmpfield(:ncol,i) = surface_state2d(i)%sols(:ncol)
      end do
      call write_field_from_chunk(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
       ncol = get_ncols_p(i)
	tmpfield(:ncol,i) = surface_state2d(i)%soll(:ncol)
      end do
      call write_field_from_chunk(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
       ncol = get_ncols_p(i)
	tmpfield(:ncol,i) = surface_state2d(i)%solsd(:ncol)
      end do
      call write_field_from_chunk(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
       ncol = get_ncols_p(i)
	tmpfield(:ncol,i) = surface_state2d(i)%solld(:ncol)
      end do
      call write_field_from_chunk(nrg,1,1,1,tmpfield)
      call write_field_from_chunk(nrg,1,1,1,trefmxav)
      call write_field_from_chunk(nrg,1,1,1,trefmnav)
      call write_field_from_chunk(nrg,1,1,1,icefrac)
      call write_field_from_chunk(nrg,1,1,1,Focn)
      call write_field_from_chunk(nrg,1,1,1,tsocn)
      call write_field_from_chunk(nrg,1,1,1,frzmlt)
#ifdef CRM
      call write_field_from_chunk(nrg,1,crm_nx*crm_ny*crm_nz,1,u_crm)
      call write_field_from_chunk(nrg,1,crm_nx*crm_ny*crm_nz,1,v_crm)
      call write_field_from_chunk(nrg,1,crm_nx*crm_ny*crm_nz,1,w_crm)
      call write_field_from_chunk(nrg,1,crm_nx*crm_ny*crm_nz,1,t_crm)
      call write_field_from_chunk(nrg,1,crm_nx*crm_ny*crm_nz,1,q_crm)
      call write_field_from_chunk(nrg,1,crm_nx*crm_ny*crm_nz,1,qn_crm)
      call write_field_from_chunk(nrg,1,crm_nx*crm_ny*crm_nz,1,qp_crm)
      call write_field_from_chunk(nrg,1,crm_nx*crm_ny*crm_nz,1,qrs_crm)
      call write_field_from_chunk(nrg,1,crm_nx*crm_ny*crm_nz,1,qrl_crm)
      call write_field_from_chunk(nrg,1,nrad_buffer,1,rad_buffer)
      call write_field_from_chunk(nrg,1,pver,1,qrs1)
      call write_field_from_chunk(nrg,1,pver,1,qrl1)
      call write_field_from_chunk(nrg,1,crm_nx*crm_ny,1,fsds_crm)
      call write_field_from_chunk(nrg,1,crm_nx*crm_ny,1,fsns_crm)
      call write_field_from_chunk(nrg,1,crm_nx*crm_ny,1,fsdsc_crm)
      call write_field_from_chunk(nrg,1,crm_nx*crm_ny,1,fsntoa_crm)
      call write_field_from_chunk(nrg,1,crm_nx*crm_ny,1,fsntoac_crm)
      call write_field_from_chunk(nrg,1,crm_nx*crm_ny,1,fsutoa_crm)
      call write_field_from_chunk(nrg,1,crm_nx*crm_ny,1,flwds_crm)
      call write_field_from_chunk(nrg,1,crm_nx*crm_ny,1,flns_crm)
      call write_field_from_chunk(nrg,1,crm_nx*crm_ny,1,flnsc_crm)
      call write_field_from_chunk(nrg,1,crm_nx*crm_ny,1,flut_crm)
      call write_field_from_chunk(nrg,1,crm_nx*crm_ny,1,flutc_crm)
#endif

#ifdef QRLDAMP
      call write_field_from_chunk (nrg, 1, pver*48,1,qrl_dailymean_buffer)
      call write_field_from_chunk (nrg,1,pver,1,qrl_dailymean)
      call write_field_from_chunk (nrg,1,pver,1,qrl_anncycle_int)
#endif
!   allocate (qrl_dailymean_buffer(pcols,pver,48,begchunk:endchunk))

#if ( defined COUP_CSM )
      call write_restart_ccsm ()
#endif



!
!-----------------------------------------------------------------------
! Write the abs/ems restart dataset if necessary    
!-----------------------------------------------------------------------
!
      if (aeres) then
         if (masterproc) then
            fname = interpret_filename_spec( rafilename_spec )
            pname = trim(get_archivedir('rest'))//fname
            call opnfil(fname, nrg2, 'u')
            write(nrg,iostat=ioerr) pname
            if (ioerr /= 0 ) then
               write (6,*) 'WRITE ioerror ',ioerr,' on i/o unit = ',nrg
               call endrun
            end if
         endif

         call write_field_from_chunk(nrg2, 1, pverp*pverp,1, abstot_3d(1,1,1,begchunk))
         call write_field_from_chunk(nrg2, 1, pver*4,     1, absnxt_3d(1,1,1,begchunk))
         call write_field_from_chunk(nrg2, 1, pverp,      1, emstot_3d(1,1,begchunk))


         if (masterproc) then
            close(nrg2)
            call putfil (fname, pname, mss_wpass, mss_irt, (.not. nlend) )
         end if
      end if


      
      return
   end subroutine write_restart_physics

!#######################################################################

   subroutine read_restart_physics (nrg, nrg2, aeres )

      use phys_buffer, only: pbuf_allocate, pbuf_read_restart
      use phys_grid,   only: get_ncols_p
!
! Arguments
!
      integer, intent(in) :: nrg
      integer, intent(in) :: nrg2

      logical, intent(in) :: aeres
!
! Local workspace
!
      real(r8) tmpfield(pcols,begchunk:endchunk)
      real(r8) tmpfield3d(pcols,plevmx,begchunk:endchunk)
      integer i                 ! loop index
      integer n3tmp             ! timestep index
      character*80  locfn       ! Local filename
      integer ioerr             ! I/O status
      integer :: ncol           ! number of columns in a chunk
!
! Buffer module variables
!
      call pbuf_allocate('global')
      call initialize_buffer ()

      call read_chunk_from_field(nrg,1,1,1,pblht)
      call read_chunk_from_field(nrg,1,1,1,tpert)
      call read_chunk_from_field(nrg,1,pver,1,qrs)
      call read_chunk_from_field(nrg,1,pver,1,qrl)
      call read_chunk_from_field(nrg,1,pcnst+pnats,1,qpert)

! Physics buffer
      call pbuf_read_restart(nrg)
!
! Comsrf module variables
!
      call initialize_comsrf
#if (! defined COUP_CSM)
      call read_chunk_from_field(nrg,1,1,1,fsnt)
#endif
      call read_chunk_from_field(nrg,1,1,1,fsns)
#if (! defined COUP_CSM)
      call read_chunk_from_field(nrg,1,1,1,flnt)
      call read_chunk_from_field(nrg,1,1,1,flns)
#endif
      call read_chunk_from_field(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
        ncol = get_ncols_p(i)
	srfflx_state2d(i)%asdir(:ncol) = tmpfield(:ncol,i)
      end do
      call read_chunk_from_field(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
        ncol = get_ncols_p(i)
	srfflx_state2d(i)%asdif(:ncol) = tmpfield(:ncol,i)
      end do
      call read_chunk_from_field(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
        ncol = get_ncols_p(i)
	srfflx_state2d(i)%aldir(:ncol) = tmpfield(:ncol,i)
      end do
      call read_chunk_from_field(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
        ncol = get_ncols_p(i)
	srfflx_state2d(i)%aldif(:ncol) = tmpfield(:ncol,i)
      end do

#if (! defined COUP_CSM)
      call read_chunk_from_field(nrg,1,1,1,asdirice)
      call read_chunk_from_field(nrg,1,1,1,asdifice)
      call read_chunk_from_field(nrg,1,1,1,aldirice)
      call read_chunk_from_field(nrg,1,1,1,aldifice)
      call read_chunk_from_field(nrg,1,1,1,asdirocn)
      call read_chunk_from_field(nrg,1,1,1,asdifocn)
      call read_chunk_from_field(nrg,1,1,1,aldirocn)
      call read_chunk_from_field(nrg,1,1,1,aldifocn)
      call read_chunk_from_field(nrg,1,1,1,asdirlnd)
      call read_chunk_from_field(nrg,1,1,1,asdiflnd)
      call read_chunk_from_field(nrg,1,1,1,aldirlnd)
      call read_chunk_from_field(nrg,1,1,1,aldiflnd)
      call read_chunk_from_field(nrg,1,1,1,lwuplnd)
      call read_chunk_from_field(nrg,1,1,1,lwupocn)
      call read_chunk_from_field(nrg,1,1,1,lwupice)
      call read_chunk_from_field(nrg,1,1,1,tsice)
#endif
      call read_chunk_from_field(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
        ncol = get_ncols_p(i)
	srfflx_state2d(i)%lwup(:ncol) = tmpfield(:ncol,i)
      end do
      call read_chunk_from_field(nrg,1,1,1,landfrac)
#ifdef COUP_SOM
      call gather_chunk_to_field(1,1,1,plon,landfrac,landfrac_field)   ! all nodes need global landfrac
      call read_chunk_from_field(nrg,1,1,1,aice)
#endif
      call read_chunk_from_field(nrg,1,1,1,ocnfrac)
      call read_chunk_from_field(nrg,1,1,1,icefrac)
      call read_chunk_from_field(nrg,1,1,1,landm)
      call read_chunk_from_field(nrg,1,1,1,sgh)
      call read_chunk_from_field(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
        ncol = get_ncols_p(i)
	srfflx_state2d(i)%ts(:ncol) = tmpfield(:ncol,i)
      end do
      call read_chunk_from_field(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
        ncol = get_ncols_p(i)
	srfflx_state2d(i)%sst(:ncol) = tmpfield(:ncol,i)
      end do
      call read_chunk_from_field(nrg,1,plevmx,1,tmpfield3d)
      do i=begchunk,endchunk
        ncol = get_ncols_p(i)
	surface_state2d(i)%tssub(:ncol,:) = tmpfield3d(:ncol,:,i)
      end do
      call read_chunk_from_field(nrg,1,1,1,sicthk)
      call read_chunk_from_field(nrg,1,1,1,snowhland)
#if (! defined COUP_CSM)
      call read_chunk_from_field(nrg,1,1,1,snowhice)
#endif
      call read_chunk_from_field(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
        ncol = get_ncols_p(i)
	surface_state2d(i)%flwds(:ncol) = tmpfield(:ncol,i)
      end do
      call read_chunk_from_field(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
        ncol = get_ncols_p(i)
	surface_state2d(i)%sols(:ncol) = tmpfield(:ncol,i)
      end do
      call read_chunk_from_field(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
        ncol = get_ncols_p(i)
	surface_state2d(i)%soll(:ncol) = tmpfield(:ncol,i)
      end do
      call read_chunk_from_field(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
        ncol = get_ncols_p(i)
	surface_state2d(i)%solsd(:ncol) = tmpfield(:ncol,i)
      end do
      call read_chunk_from_field(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
        ncol = get_ncols_p(i)
	surface_state2d(i)%solld(:ncol) = tmpfield(:ncol,i)
      end do
      call read_chunk_from_field(nrg,1,1,1,trefmxav)
      call read_chunk_from_field(nrg,1,1,1,trefmnav)
      call read_chunk_from_field(nrg,1,1,1,icefrac)
      call read_chunk_from_field(nrg,1,1,1,Focn)
      call read_chunk_from_field(nrg,1,1,1,tsocn)
      call read_chunk_from_field(nrg,1,1,1,frzmlt)
#ifdef CRM
      call read_chunk_from_field(nrg,1,crm_nx*crm_ny*crm_nz,1,u_crm)
      call read_chunk_from_field(nrg,1,crm_nx*crm_ny*crm_nz,1,v_crm)
      call read_chunk_from_field(nrg,1,crm_nx*crm_ny*crm_nz,1,w_crm)
      call read_chunk_from_field(nrg,1,crm_nx*crm_ny*crm_nz,1,t_crm)
      call read_chunk_from_field(nrg,1,crm_nx*crm_ny*crm_nz,1,q_crm)
      call read_chunk_from_field(nrg,1,crm_nx*crm_ny*crm_nz,1,qn_crm)
      call read_chunk_from_field(nrg,1,crm_nx*crm_ny*crm_nz,1,qp_crm)
      call read_chunk_from_field(nrg,1,crm_nx*crm_ny*crm_nz,1,qrs_crm)
      call read_chunk_from_field(nrg,1,crm_nx*crm_ny*crm_nz,1,qrl_crm)
      call read_chunk_from_field(nrg,1,nrad_buffer,1,rad_buffer)
      call read_chunk_from_field(nrg,1,pver,1,qrs1)
      call read_chunk_from_field(nrg,1,pver,1,qrl1)
      call read_chunk_from_field(nrg,1,crm_nx*crm_ny,1,fsds_crm)
      call read_chunk_from_field(nrg,1,crm_nx*crm_ny,1,fsns_crm)
      call read_chunk_from_field(nrg,1,crm_nx*crm_ny,1,fsdsc_crm)
      call read_chunk_from_field(nrg,1,crm_nx*crm_ny,1,fsntoa_crm)
      call read_chunk_from_field(nrg,1,crm_nx*crm_ny,1,fsntoac_crm)
      call read_chunk_from_field(nrg,1,crm_nx*crm_ny,1,fsutoa_crm)
      call read_chunk_from_field(nrg,1,crm_nx*crm_ny,1,flwds_crm)
      call read_chunk_from_field(nrg,1,crm_nx*crm_ny,1,flns_crm)
      call read_chunk_from_field(nrg,1,crm_nx*crm_ny,1,flnsc_crm)
      call read_chunk_from_field(nrg,1,crm_nx*crm_ny,1,flut_crm)
      call read_chunk_from_field(nrg,1,crm_nx*crm_ny,1,flutc_crm)
#endif

#ifdef QRLDAMP
      call read_chunk_from_field (nrg,1,pver*48,1,qrl_dailymean_buffer)
      call read_chunk_from_field (nrg,1,pver,1,qrl_dailymean)
      call read_chunk_from_field (nrg,1,pver,1,qrl_anncycle_int)

#endif


#if ( defined COUP_CSM )
      call initialize_ccsm_msg ()
      call read_restart_ccsm ()
#endif
!
!-----------------------------------------------------------------------
! Read the abs/ems restart dataset if necessary    
!-----------------------------------------------------------------------
!
      call initialize_radbuffer ()
      if (aeres) then
         if (masterproc) then
            read(nrg,iostat=ioerr) pname
            if (ioerr /= 0 ) then
               write (6,*) 'READ ioerror ',ioerr,' on i/o unit = ',nrg
               call endrun
            end if
            call getfil (pname, locfn)
            call opnfil (locfn, nrg2, 'u')
         endif

         call read_chunk_from_field(nrg2, 1, pverp*pverp,1,abstot_3d(1,1,1,begchunk))
         call read_chunk_from_field(nrg2, 1, pver*4,     1,absnxt_3d(1,1,1,begchunk))
         call read_chunk_from_field(nrg2, 1, pverp,      1,emstot_3d(1,1,begchunk))

         if (masterproc) close(nrg2)
      end if
      
      return
   end subroutine read_restart_physics

   character(len=256) function get_abs_restart_filepath ( )
!
! Return the full filepath to the abs-ems restart file
!
     get_abs_restart_filepath = pname
   end function get_abs_restart_filepath

end module restart_physics
