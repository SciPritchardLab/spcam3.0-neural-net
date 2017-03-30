#include <params.h>
#ifdef CRM
module crmics
   use crmdims, only: crm_nx, crm_ny, crm_nz, crm_dx, crm_dy
   use runtime_opts, only: crminitread,crmsavechunks
   use buffer
   use phys_grid
   use pmgrid, only: plev

implicit none


contains

   subroutine read_crm_ics ()
 
   implicit none   

#include <comlun.h>

   integer :: m
   integer :: ucrmid(crmsavechunks)
   integer :: vcrmid(crmsavechunks)
   integer :: tcrmid(crmsavechunks)
   integer :: wcrmid(crmsavechunks)
   integer :: qcrmid(crmsavechunks)
   integer :: qpcrmid(crmsavechunks)
   integer :: qncrmid(crmsavechunks)
   integer :: qrscrmid(crmsavechunks)
   integer :: qrlcrmid(crmsavechunks)
   character(1) numb

      integer :: fsdscrmid,fsnscrmid,fsdsccrmid
      integer :: fsntoacrmid,fsntoaccrmid
      integer :: fsutoacrmid,flwdscrmid,flnscrmid
      integer :: flutccrmid,flnsccrmid,flutcrmid
      integer :: qrs1id,qrl1id,crmradbfrid
  
  integer start3d(4)             ! starting index array for nf_put_vara (3d vars)
  integer count3d(4)             ! count array for nf_put_vara (3d vars)
     integer start3dcrmradbfr(4)
      integer count3dcrmradbfr(4)
   
    real(r8), allocatable :: arr3d_lonlat(:,:,:) 
   real(r8), allocatable :: arr3d_colchnk(:,:,:) 
    real(r8), allocatable :: arr3dcrmradbfr_lonlat(:,:,:) 
   real(r8), allocatable :: arr3dcrmradbfr_colchnk(:,:,:) 

   start3d(1) = 1
   start3d(2) = 1
   start3d(3) = 1
   start3d(4) = 1
   start3dcrmradbfr = start3d
   count3d(1) = plon
   count3d(2) = plev
   count3d(3) = plat
    count3d(4) = 1
   count3dcrmradbfr = count3d
   count3dcrmradbfr(2) = nrad_buffer


    if (crminitread ) then
       do m=1,crmsavechunks
           if (masterproc) then
             write (6,*) 'XYXY in crminitread ncid_ini=',ncid_ini
             write(numb,'(i1)') m
             numb = adjustl(numb)
            call wrap_inq_varid (ncid_ini,'U_CRM'//trim(numb),ucrmid(m))
            call wrap_inq_varid (ncid_ini,'V_CRM'//trim(numb),vcrmid(m))
            call wrap_inq_varid (ncid_ini,'T_CRM'//trim(numb),tcrmid(m))
            call wrap_inq_varid (ncid_ini,'W_CRM'//trim(numb),wcrmid(m))
            call wrap_inq_varid (ncid_ini,'Q_CRM'//trim(numb),qcrmid(m))
            call wrap_inq_varid (ncid_ini,'QP_CRM'//trim(numb),qpcrmid(m))
            call wrap_inq_varid (ncid_ini,'QN_CRM'//trim(numb),qncrmid(m))
            call wrap_inq_varid (ncid_ini,'QRS_CRM'//trim(numb),qrscrmid(m))
            call wrap_inq_varid (ncid_ini,'QRL_CRM'//trim(numb),qrlcrmid(m))
           end if
       end do
      call reconstitute_crm3d_icvar (ucrmid(:), u_crm ) !
       call reconstitute_crm3d_icvar (vcrmid(:), v_crm ) !
       call reconstitute_crm3d_icvar (tcrmid(:), t_crm ) !
       call reconstitute_crm3d_icvar (wcrmid(:), w_crm ) !
       call reconstitute_crm3d_icvar (qcrmid(:), q_crm ) !
       call reconstitute_crm3d_icvar (qncrmid(:), qn_crm ) !
       call reconstitute_crm3d_icvar (qpcrmid(:), qp_crm ) !
       call reconstitute_crm3d_icvar (qrscrmid(:), qrs_crm ) !
       call reconstitute_crm3d_icvar (qrlcrmid(:), qrl_crm ) !

       ! Now do the 2D variables:

     if (masterproc) then
       call wrap_inq_varid (ncid_ini,'FSDS_CRM',fsdscrmid)
       call wrap_inq_varid (ncid_ini,'FSNS_CRM',fsnscrmid)
       call wrap_inq_varid (ncid_ini,'FSDSC_CRM',fsdsccrmid)
       call wrap_inq_varid (ncid_ini,'FSNTOA_CRM',fsntoacrmid)
       call wrap_inq_varid (ncid_ini,'FSNTOAC_CRM',fsntoaccrmid)
       call wrap_inq_varid (ncid_ini,'FSUTOA_CRM',fsutoacrmid)
       call wrap_inq_varid (ncid_ini,'FLWDS_CRM',flwdscrmid)
       call wrap_inq_varid (ncid_ini,'FLNS_CRM',flnscrmid)
       call wrap_inq_varid (ncid_ini,'FLNSC_CRM',flnsccrmid)
       call wrap_inq_varid (ncid_ini,'FLUT_CRM',flutcrmid)
       call wrap_inq_varid (ncid_ini,'FLUTC_CRM',flutccrmid)
     end if

     call reconstitute_crm2d_icvar (fsdscrmid, fsds_crm)
     call reconstitute_crm2d_icvar (fsnscrmid, fsns_crm)
     call reconstitute_crm2d_icvar (fsdsccrmid, fsdsc_crm)
     call reconstitute_crm2d_icvar (fsntoacrmid, fsntoa_crm)
     call reconstitute_crm2d_icvar (fsntoaccrmid, fsntoac_crm)
     call reconstitute_crm2d_icvar (fsutoacrmid, fsutoa_crm)
     call reconstitute_crm2d_icvar (flwdscrmid, flwds_crm)
     call reconstitute_crm2d_icvar (flnscrmid, flns_crm)
     call reconstitute_crm2d_icvar (flnsccrmid, flnsc_crm)
     call reconstitute_crm2d_icvar (flutcrmid, flut_crm)
     call reconstitute_crm2d_icvar (flutccrmid, flutc_crm)

     ! Now do the variables without CRM interior dimensions to unfold:
      if (masterproc) then
        allocate (arr3d_lonlat (plon,plev,plat) ) ! This seems consistent with history.F90
        call wrap_inq_varid (ncid_ini,'QRS1',qrs1id)
        write (6,*) 'MCHECK1 read var', qrs1id
        call wrap_get_vara_realx(ncid_ini,qrs1id,start3d,count3d,arr3d_lonlat)
      end if  
     call scatter_field_to_chunk(1,plev,1,plon,arr3d_lonlat,qrs1)
     if (masterproc) then
        call wrap_inq_varid (ncid_ini,'QRL1',qrl1id)
        write (6,*) 'MCHECK2 read var', qrl1id
        call wrap_get_vara_realx(ncid_ini,qrl1id,start3d,count3d,arr3d_lonlat)
     end if
     call scatter_field_to_chunk(1,plev,1,plon,arr3d_lonlat,qrl1)
     
     if (masterproc) then
       deallocate(arr3d_lonlat)
       allocate (arr3dcrmradbfr_lonlat (plon,nrad_buffer,plat))
       call wrap_inq_varid (ncid_ini,'CRMRADBFR',crmradbfrid)
        write (6,*) 'MCHECK3 read var', crmradbfrid
       call wrap_get_vara_realx(ncid_ini,crmradbfrid,start3dcrmradbfr,count3dcrmradbfr,arr3dcrmradbfr_lonlat)
      endif 
      call scatter_field_to_chunk(1,nrad_buffer,1,plon,arr3dcrmradbfr_lonlat,rad_buffer)
      if (masterproc) then
        deallocate(arr3dcrmradbfr_lonlat)
      end if

!       if (masterproc) then
!          write (6,*) 'MCHECK INIT: Q_CRM(3,22,1,:,begchunk) = ', q_crm(3,22,1,:,begchunk)
!          write (6,*) 'MCHECK INIT: Q_CRM(7,14,1,:,endchunk) = ', q_crm(7,14,1,:,endchunk)
!       end if
!      MSP NOTE this sanity check passed: matched the output in history.F90
       end if
end subroutine read_crm_ics

subroutine reconstitute_crm3d_icvar(crmvarid,crmvar_reconstituted)

   use runtime_opts, only: crminitread,crmsavechunks
   use phys_grid

   implicit none

#include <comlun.h>

   integer, intent(in) :: crmvarid(crmsavechunks)
   real(r8), intent(inout) :: crmvar_reconstituted(pcols,crm_nx,crm_ny,crm_nz,begchunk:endchunk)
   integer start3dcrm(4)             ! starting index array for nf_put_vara
   integer count3dcrm(4)
   data start3dcrm/1,1,1,1/
   integer :: chunkdiml,m
   real(r8), allocatable :: arr3dcrmchnk_lonlat(:,:,:) ! For pulling individual chunked variables
   real(r8), allocatable :: arr3dcrmchnk_colchnk(:,:,:) ! For individual chunks remapped from lon/lat to physics grid domain decomp (col/chnk)
   integer :: k1,k2,icol,cc,icrmx,icrmy,kk
   character(1) numb

   chunkdiml = crm_nx*crm_ny*crm_nz/crmsavechunks
   count3dcrm(1) = plon
   count3dcrm(2) = chunkdiml
   count3dcrm(3) = plat
   count3dcrm(4) = 1
   allocate (arr3dcrmchnk_lonlat (plon,chunkdiml,plat) ) ! This seems consistent with history.F90
   allocate (arr3dcrmchnk_colchnk (pcols,chunkdiml,begchunk:endchunk) )

   do m=1,crmsavechunks
   if (masterproc) then ! Only masterproc reads the IC file; other procs see scatter command though and receive the results.  
      write (numb,'(i1)') m
      call wrap_get_vara_realx(ncid_ini, crmvarid(m), start3dcrm,count3dcrm,arr3dcrmchnk_lonlat)
   end if
   ! Next scatter from lon,lat --> column, chunk of physics space.
   ! Recall the reason we collapsed the interior CRM dimensions was to leverage this function which does that for us (but only for a 3D variable):
   call scatter_field_to_chunk(1,chunkdiml,1,plon,arr3dcrmchnk_lonlat,arr3dcrmchnk_colchnk)
   ! And I suspect the reason we "chunked" is due to finite buffer size for above. 
   ! Now put the chunk where it belongs and unfold the CRM dimensions:
   k1 = 1+(m-1)*crm_nz/crmsavechunks
   k2 = m*crm_nz/crmsavechunks
   do icol=1,pcols
     cc=0
     do icrmx=1,crm_nx
       do icrmy=1,crm_ny
         do kk=k1,k2
           cc=cc+1
           crmvar_reconstituted(icol,icrmx,icrmy,kk,begchunk:endchunk) = arr3dcrmchnk_colchnk (icol,cc,begchunk:endchunk)
          ! The above expression is consistent to what was in history.F90.
         end do
       end do
     end do
   end do
   end do ! over chunks.
   deallocate (arr3dcrmchnk_lonlat)
   deallocate (arr3dcrmchnk_colchnk)
 end subroutine reconstitute_crm3d_icvar

subroutine reconstitute_crm2d_icvar (crm2dvarid,crm2dvar_reconstituted)
                                                       
   use runtime_opts, only: crminitread
   use phys_grid

   implicit none

#include <comlun.h>

   integer, intent(in) :: crm2dvarid
   real(r8), intent(inout) :: crm2dvar_reconstituted(pcols,crm_nx,crm_ny,begchunk:endchunk)
   integer start3dcrm(4)             ! starting index array for nf_put_vara
   integer count3dcrm(4)
   data start3dcrm/1,1,1,1/
   integer :: chunkdiml,m

   real(r8), allocatable :: arr3dcrmchnk_lonlat(:,:,:) ! For pulling individual chunked variables
   real(r8), allocatable :: arr3dcrmchnk_colchnk(:,:,:) ! For individual chunks remapped from lon/lat to physics grid domain decomp (col/chnk)
   integer :: icol,cc,icrmx,icrmy,kk
   character(1) numb

   chunkdiml = crm_nx*crm_ny
   count3dcrm(1) = plon
   count3dcrm(2) = chunkdiml
   count3dcrm(3) = plat
   count3dcrm(4) = 1
   allocate (arr3dcrmchnk_lonlat (plon,chunkdiml,plat) ) ! This seems consistent with history.F90
   allocate (arr3dcrmchnk_colchnk (pcols,chunkdiml,begchunk:endchunk) )

   if (masterproc) then ! Only masterproc reads the IC file; other procs see scatter command though and receive the results.  
      call wrap_get_vara_realx(ncid_ini, crm2dvarid, start3dcrm,count3dcrm,arr3dcrmchnk_lonlat)
   end if
   ! Next scatter from lon,lat --> column, chunk of physics space.
   ! Recall the reason we collapsed the interior CRM dimensions was to leverage this function which does that for us (but only for a 3D variable):
   call scatter_field_to_chunk(1,chunkdiml,1,plon,arr3dcrmchnk_lonlat,arr3dcrmchnk_colchnk)
   ! And I suspect the reason we "chunked" is due to finite buffer size for above. 
   ! Now put the chunk where it belongs and unfold the CRM dimensions:
   do icol=1,pcols
     cc=0
     do icrmx=1,crm_nx
       do icrmy=1,crm_ny
          cc=cc+1
          crm2dvar_reconstituted(icol,icrmx,icrmy,begchunk:endchunk) = arr3dcrmchnk_colchnk (icol,cc,begchunk:endchunk)
         ! The above expression is consistent to what was in history.F90.
       end do
     end do
   end do
   deallocate (arr3dcrmchnk_lonlat)
   deallocate (arr3dcrmchnk_colchnk)
 end subroutine reconstitute_crm2d_icvar

end module crmics
#endif
