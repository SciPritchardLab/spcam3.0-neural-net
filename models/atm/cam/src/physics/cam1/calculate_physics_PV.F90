#ifdef PVBUDGET
subroutine calculate_physics_PV (state,pv1,pv2,pv3)

!----------------------------------------------------------------------- 
! 
! Author: M. S. Pritchard
!
! A port of an algorithm contributed to NCL by D. Shea
! 
! Purpose: 
! Estimate the three terms in this equation for PV:
!   pv    = -G*(vr*dthdp - (dthdx*dvdp-dthdy*dudp) )
! Units: N/kg * (K/m*m/s/Pa) = N/kg*K/(N/m2) = (K m2/kg)/s

! WARNING: assumes that state%t,u,v,pmid are all up to date
! (when running in the SLD dycore, it's possible t will not be depending on where in the code we are applying this)
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use physics_types,      only: physics_state
   use physconst,       only: cappa,gravit
   use phys_grid, only: get_ncols_p,clon_p, clat_p, scatter_field_to_chunk, gather_chunk_to_field ! Full physics grid lons & lats (radians)
   use pmgrid, only: beglat, endlat, plon, plat, masterproc
   use ppgrid, only: pcols, pver, begchunk, endchunk

   implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments
  type(physics_state), intent(inout), dimension (begchunk:endchunk) :: state
! Output arguments
  real(r8), intent(out) :: pv1(pcols,pver,begchunk:endchunk)
  real(r8), intent(out) :: pv2(pcols,pver,begchunk:endchunk)
  real(r8), intent(out) :: pv3(pcols,pver,begchunk:endchunk)

! Local variables:
 
  real (r8) :: theta(pcols,pver)
  real (r8) :: dthetadp(pcols,pver),dudp(pcols,pver),dvdp(pcols,pver)

  real (r8) :: dthetadx(pcols,pver,begchunk:endchunk), dthetady(pcols,pver,begchunk:endchunk), relvort(pcols,pver,begchunk:endchunk)

  integer :: ncol,c,i,j,k
  real (r8) :: re,twopi
  real(r8), allocatable :: arr3d(:,:,:) 
  real(r8), allocatable ::  arr3dchnked (:,:,:) 

  real(r8), allocatable :: tglobal(:,:,:) 
  real(r8), allocatable :: uglobal(:,:,:) 
  real(r8), allocatable :: vglobal(:,:,:) 
  real(r8), allocatable :: thetaglobal(:,:,:) 
  real(r8), allocatable :: pmidglobal(:,:,:) 
  real(r8), allocatable :: dthetadxglobal(:,:,:) 
  real(r8), allocatable :: dthetadyglobal(:,:,:) 
  real(r8), allocatable :: relvortglobal(:,:,:) 



  real (r8) :: dx, dy, dx_radians,dy_radians,duy,dvx,dthetax,dthetay,dvdx,dudy

   if (masterproc) then
         allocate (arr3d(plon,pver,plat))
         allocate (tglobal(plon,pver,plat))
         allocate (thetaglobal(plon,pver,plat))
         allocate (uglobal(plon,pver,plat))
         allocate (vglobal(plon,pver,plat))
       allocate (pmidglobal(plon,pver,plat))
         allocate (dthetadxglobal(plon,pver,plat))
         allocate (dthetadyglobal(plon,pver,plat))
         allocate (relvortglobal(plon,pver,plat))
   end if
!         allocate (arr3d(plon,pver,beglat:endlat))
!         allocate (tglobal(plon,pver,beglat:endlat))
!         allocate (uglobal(plon,pver,beglat:endlat))
!         allocate (vglobal(plon,pver,beglat:endlat))
!         allocate (dthetadxglobal(plon,pver,beglat:endlat))
!         allocate (dthetadyglobal(plon,pver,beglat:endlat))
!         allocate (relvortglobal(plon,pver,beglat:endlat))
!   end if

   allocate (arr3dchnked (pcols,pver,begchunk:endchunk)) 


   twopi = 2.*3.141592654

   ! Master processor does all the horizontal derivatives (necessary):
 
   ! Gather global fields of U on master processor into (uglobal,vglobal,tglobal,pmidglobal) 
     do c=begchunk,endchunk
       arr3dchnked (1:pcols,1:pver,c) = state(c)%u(1:pcols,1:pver)
     end do
     call gather_chunk_to_field(1,pver,1,plon,arr3dchnked(1:pcols,1:pver,begchunk:endchunk),uglobal)

   ! and V
     do c=begchunk,endchunk
       arr3dchnked (1:pcols,1:pver,c) = state(c)%v(1:pcols,1:pver)
     end do
     call gather_chunk_to_field(1,pver,1,plon,arr3dchnked(1:pcols,1:pver,begchunk:endchunk),vglobal)
   ! and T  
     do c=begchunk,endchunk
       arr3dchnked (1:pcols,1:pver,c) = state(c)%t(1:pcols,1:pver)
     end do
     call gather_chunk_to_field(1,pver,1,plon,arr3dchnked(1:pcols,1:pver,begchunk:endchunk),tglobal)

   ! and p (to convert T to theta)  
     do c=begchunk,endchunk
       arr3dchnked (1:pcols,1:pver,c) = state(c)%pmid(1:pcols,1:pver)
     end do
     call gather_chunk_to_field(1,pver,1,plon,arr3dchnked(1:pcols,1:pver,begchunk:endchunk),pmidglobal)

   ! Master processor computes the horizontal derivatives 

     if (masterproc) then 
       ! Calculate theta (assumes t is up to date in state, need to be careful to check this!).
       do i=1,plon
        do j=1,plat
          do k=1,pver
            thetaglobal(i,k,j) = tglobal (i,k,j)*(100000./pmidglobal(i,k,j))**cappa
          end do
        end do
       end do 
  

 
       re = 6.37122e06 ! Mean earth radius (m) from cam_diagnostics.ncl
       do i=1,plon
         do j=1,plat
           if (i .eq. 1) then
             dx_radians = clon_p(2,j)-clon_p(plon,j)
          else if (i .eq. plon) then
             dx_radians = clon_p(1,j) - clon_p(plon-1,j)
           else
             dx_radians = clon_p(i+1,j)-clon_p(i-1,j)
           end if
           if (dx_radians .ge. twopi) then
             dx_radians = dx_radians - twopi
           else if (dx_radians .le. -twopi) then
             dx_radians = dx_radians + twopi
           end if

           do k=1,pver 
             theta = tglobal (i,k,j)*(100000./pmidglobal(i,k,j))**cappa

              if (j .eq. 1) then
                dy_radians = clat_p(2)-clat_p(1)
                duy = uglobal(i,k,2) - uglobal (i,k,1)
                dthetay = thetaglobal(i,k,2) - thetaglobal (i,k,1)
              else if (j .eq. plat) then
                dy_radians = clat_p(plat) - clat_p(plat-1)
                duy = uglobal(i,k,plat) - uglobal(i,k,plat-1)
                dthetay = thetaglobal(i,k,plat) - thetaglobal(i,k,plat-1)
              else
                dy_radians = clat_p(j+1)-clat_p(j-1)
                duy = uglobal(i,k,j+1) - uglobal (i,k,j-1)
                dthetay = thetaglobal(i,k,j+1) - thetaglobal (i,k,j-1)
              end if

              dx = re*cos(dy_radians)*dx_radians
              dy = re*dy_radians
  
              dudy = duy/dy
              dthetadyglobal(i,k,j) = dthetay/dy

              if (i .eq. 1) then
                dvx = uglobal(2,k,j) - vglobal(plon,k,j)
                dthetax = uglobal(2,k,j) - thetaglobal(plon,k,j)
              else if (i.eq. plon) then
                dvx = vglobal(1,k,j) - vglobal(plon-1,k,j)
                dthetax = thetaglobal(1,k,j) - thetaglobal(plon-1,k,j)
              else
                dvx = vglobal(i+1,k,j) - vglobal(i-1,k,j)
                dthetax = thetaglobal(i+1,k,j) - thetaglobal(i-1,k,j)
              end if
              
              dvdx = dvx/dx

              dthetadxglobal(i,k,j) = dthetax/dx
              relvortglobal(i,k,j) = dvdx - dudy
              end do ! k, level 
           end do ! j, latitude
        end do ! i, longitude
     end if  ! if masterproc
 
     ! Masterproc scatters horizontal derivatives to physical domain decomp "chunks":

     call scatter_field_to_chunk(1,pver,1,plon,dthetadxglobal,dthetadx(:pcols,:pver,begchunk:endchunk))
     call scatter_field_to_chunk(1,pver,1,plon,dthetadyglobal,dthetady(:pcols,:pver,begchunk:endchunk))
     call scatter_field_to_chunk(1,pver,1,plon,relvortglobal,relvort(:pcols,:pver,begchunk:endchunk))

! Below here, we are back to working in local physics chunk portion:
    do c = begchunk,endchunk
      ncol = get_ncols_p (c)
      do i=1,ncol
        do k=1,pver
          theta (i,k) = state(c)%t(i,k)*(100000./state(c)%pmid(i,k))**cappa
        end do
      end do

    ! Use central finite differencing to estimate the pressure gradient terms
    ! dthdp, dudp, dvdp
      do i=1,ncol
        dthetadp(i,1) = (theta(i,2)-theta(i,1))/(state(c)%pmid(i,2) - state(c)%pmid(i,1))
        dthetadp(i,pver) = (theta(i,pver)-theta(i,pver-1))/(state(c)%pmid(i,pver) - state(c)%pmid(i,pver-1))
        dudp(i,1) = (state(c)%u(i,2)-state(c)%u(i,1))/(state(c)%pmid(i,2) - state(c)%pmid(i,1))
        dudp(i,pver) = (state(c)%u(i,pver)-state(c)%u(i,pver-1))/(state(c)%pmid(i,pver) - state(c)%pmid(i,pver-1))
        dvdp(i,1) = (state(c)%v(i,2)-state(c)%v(i,1))/(state(c)%pmid(i,2) - state(c)%pmid(i,1))
        dvdp(i,pver) = (state(c)%v(i,pver)-state(c)%v(i,pver-1))/(state(c)%pmid(i,pver) - state(c)%pmid(i,pver-1))
        do k=2,pver-1
            dthetadp(i,k) = (theta(i,k+1) - theta(i,k-1))/(state(c)%pmid(i,k+1) - state(c)%pmid(i,k-1))
            dudp(i,k) = (state(c)%u(i,k+1) - state(c)%u(i,k-1))/(state(c)%pmid(i,k+1) - state(c)%pmid(i,k-1))
            dvdp(i,k) = (state(c)%v(i,k+1) - state(c)%v(i,k-1))/(state(c)%pmid(i,k+1) - state(c)%pmid(i,k-1))
        end do
       ! Put the horizontal and vertical derivative terms together in the form of PV:
!  pv    = -G*(vr*dthdp - (dthdx*dvdp-dthdy*dudp) )
        do k=1,pver
          pv1(i,k,c) = -gravit*(relvort(i,k,c)*dthetadp(i,k))
          pv2(i,k,c) = +gravit*dthetadx(i,k,c)*dvdp(i,k)
          pv3(i,k,c) = -gravit*dthetady(i,k,c)*dudp(i,k)
          state(c)%pv(i,k) = pv1(i,k,c) + pv2(i,k,c) + pv3(i,k,c)
          if (state(c)%pv(i,k) .ne. state(c)%pv(i,k)) then
            write (6,*) 'MDEBUG NaN trap found invalid pv value at (i,kc) = ',i,k,c
            write (6,*) 'MDEBUG revort = ',relvort(i,k,c), ', dthetadp=',dthetadp(i,k)
            write (6,*) 'MDEBUG dthetadx = ',dthetadx(i,k,c), ', dthetady=',dthetady(i,k,c)
            write (6,*) 'MDEBUG dudp = ',dudp(i,k), ', dvdp=',dvdp(i,k)
            write (6,*) 'pv1=',pv1(i,k,c), ', pv2=',pv2(i,k,c),', pv3=',pv3(i,k,c)
            stop 
          end if
         end do
      end do ! icol
    end do ! c

  deallocate (arr3dchnked)

  if (masterproc) then 
    deallocate (arr3d)
    deallocate (tglobal)
    deallocate(uglobal)
    deallocate(vglobal)
    deallocate(thetaglobal)
    deallocate(pmidglobal)
    deallocate(dthetadxglobal)
    deallocate(dthetadyglobal)
    deallocate(relvortglobal)
  end if

end subroutine calculate_physics_PV
#endif
