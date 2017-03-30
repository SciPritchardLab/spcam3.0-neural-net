#include <misc.h>
#include <params.h>

subroutine sltini(dlam,    sinlam,  coslam,  uxl,     uxr, &
                  vxl,     vxr,     qxl,     qxr,     u3,  &
                  v3,      qminus,  n3m1)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Prepare the extended arrays for use in the SLT routines
!
!   1)  Fill latitude extensions.
!   2)  Fill longitude extensions.
!   3)  Compute x-derivatives
! 
! Method: 
! Computational note: The latitude loop in this routine is multitasked
! 
! Author: 
! Original version:  J. Olson
! Standardized:      J. Rosinski, June 1992
! Reviewed:          D. Williamson, P. Rasch, August 1992
! Reviewed:          D. Williamson, P. Rasch, March 1996
!
!-----------------------------------------------------------------------
!
! $Id: sltini.F90,v 1.10.4.2 2002/06/15 13:47:53 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use constituents, only: pcnst
   use rgrid
   use prognostics,  only: ptimelevels
#ifdef SPMD
   use spmd_dyn, only: cutex
#endif
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
#include <parslt.h>
!---------------------------Local parameters----------------------------
!
   integer puvpts            ! number of u/v pts in lat slice
   integer pqpts             ! number of constituent pts in lat slice
!
   parameter(puvpts = plond*plev, pqpts  = plond*plev*pcnst) 
!-----------------------------------------------------------------------
!
! Input arguments
!
   real(r8), intent(in) :: dlam(platd)          ! increment in x-direction
   real(r8), intent(in) :: sinlam(plond,platd)  ! sin(lamda)
   real(r8), intent(in) :: coslam(plond,platd)  ! cos(lamda)
   real(r8), intent(inout) :: uxl (plond,plev,      beglatex:endlatex) 
   real(r8), intent(inout) :: uxr (plond,plev,      beglatex:endlatex)  
   real(r8), intent(inout) :: vxl (plond,plev,      beglatex:endlatex) 
   real(r8), intent(inout) :: vxr (plond,plev,      beglatex:endlatex)  
   real(r8), intent(inout) :: qxl (plond,plev,pcnst,beglatex:endlatex)  
   real(r8), intent(inout) :: qxr (plond,plev,pcnst,beglatex:endlatex)

   real(r8), intent(inout) :: u3(plond, plev, beglatex:endlatex,ptimelevels)     ! u-wind component
   real(r8), intent(inout) :: v3(plond, plev, beglatex:endlatex,ptimelevels)     ! v-wind component
   real(r8), intent(inout) :: qminus(plond, plev, pcnst, beglatex:endlatex) ! moisture

   integer, intent(in) ::  n3m1                      ! time indicies
!                                    
!
!-----------------------------------------------------------------------
!
!  dlam    Length of increment in longitude grid.
!  sinlam  Sin of longitudes in global grid (model grid pts only).
!  coslam  Cos of longitudes in global grid (model grid pts only).
!  uxl     x-derivatives of u at the left  (west) edge of given interval
!  vxl     x-derivatives of v at the left  (west) edge of given interval
!  uxr     x-derivatives of u at the right (east) edge of given interval
!  vxr     x-derivatives of v at the right (east) edge of given interval
!  qxl     x-derivatives of scalar species at the left  (west) edge
!          of given interval
!  qxr     x-derivatives of scalar species at the right (east) edge
!          of given interval
!
!---------------------------Local variables-----------------------------
!
   integer m,j,k             ! index
   integer nlond
!
!------------------------------Externals--------------------------------
!
   external cubxdr,extx,extys,extyv,limdx
!
!-----------------------------------------------------------------------
!
! Fill latitude extensions beyond the southern- and northern-most
! latitudes in the global grid
!
   call t_startf ('slt_single')
   call extyv(1, plev, coslam, sinlam, u3(1,1,beglatex,n3m1), &
                                       v3(1,1,beglatex,n3m1))
   call extys(pcnst, plev    ,qminus, pcnst)
!
! Fill longitude extensions
!
   call extx(1 ,plev    ,u3(1,1,beglatex,n3m1), 1)
   call extx(1 ,plev    ,v3(1,1,beglatex,n3m1), 1)
   call extx(pcnst, plev    ,qminus, pcnst)
   call t_stopf ('slt_single')
!
! Compute x-derivatives.
!
#if ( defined SPMD )
   do j=cutex(1,iam)+numbnd,cutex(2,iam)+numbnd
#else

! The following autotasking preprocessor directives declare
! that all iterations of the loop may be done in parallel. 
! A synchronization point is introduced at the end of the loop. 
! Only variable "j" needs separate storage for each processor.

!$OMP  PARALLEL DO PRIVATE (J, K, M, NLOND)

      do j = 1,platd
#endif
         nlond = 1 + 2*nxpt + nlonex(j)
         do k=1,plev
            call cubxdr (nlond, 2, nlond-3, dlam(j), u3(1,k,j,n3m1), &
                         uxl(1,k,j), uxr(1,k,j))
            call cubxdr (nlond, 2, nlond-3, dlam(j), v3(1,k,j,n3m1), &
                         vxl(1,k,j), vxr(1,k,j))
            do m=1,pcnst
               call cubxdr (nlond, 2, nlond-3, dlam(j), qminus(1,k,m,j), &
                            qxl(1,k,m,j), qxr(1,k,m,j))
               if( plimdr )then
                  call limdx (nlond, 2, nlond-3, dlam(j), qminus(1,k,m,j), &
                              qxl(1,k,m,j), qxr(1,k,m,j))
               end if
            end do
         end do
      end do

      return
   end subroutine sltini
