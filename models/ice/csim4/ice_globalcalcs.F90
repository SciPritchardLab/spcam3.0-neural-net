#include <misc.h>

module ice_globalcalcs
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid,       only: pcols, begchunk, endchunk
   use rgrid,        only: nlon
   use commap,       only: w
   use phys_grid,    only: gather_chunk_to_field
   use pmgrid,       only: plon, plat, masterproc
   use comsrf,       only: surface_state, srfflx_parm, landfrac, aice, &
                           Focn, frzmlt, landfrac_field, sicthk, snowhice
   use ice_constants,only: ni, saltz, rLfs, Tffresh, rhofresh, Lfus, rhos, emissivity_ice, Lvap
   use phys_grid,    only: get_ncols_p, get_rlat_p, get_rlon_p
   use history,      only: outfld
#if ( defined SPMD )
   use mpishorthand
#endif

   implicit none

   CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(r8) function sea_ice_energy (srf_state, dtime, inout, deltae, deltaaice, snowfall)
      use ice_dh, only: energ
!
! Arguments
!
      type(surface_state), intent(in) :: srf_state(begchunk:endchunk)
      integer, intent(in) :: dtime
      integer, intent(in) :: inout
      real(r8), intent(inout) :: deltae(pcols,begchunk:endchunk)
      real(r8), intent(inout) :: deltaaice(pcols,begchunk:endchunk)
      real(r8), intent(in) :: snowfall(pcols,begchunk:endchunk)
!
! Local workspace
!
      real(r8) :: totenerg(pcols,begchunk:endchunk)
      real(r8) :: totenerg_field(plon,plat)
      real(r8) :: qi
      real(r8) :: wt
      real(r8) :: mean
      real(r8) :: tmp
      real(r8) :: wght
      real(r8) :: erate(pcols,begchunk:endchunk)
      real(r8) :: hs

      integer :: i, j, k, c, ncol

      character(len=8) :: varname

      call t_startf ('sea_ice_energy')
      if (inout == 1) then
         varname = 'EICEIN'
      else
         varname = 'EICEOUT'
      end if
!$OMP PARALLEL DO PRIVATE (C, I, K, NCOL, QI, HS)
      do c=begchunk,endchunk
         ncol = get_ncols_p(c)
         do i=1,ncol
            totenerg(i,c) = 0.
            if (aice(i,c) > 0.) then
               qi = 0.
               do k=1,ni
                  qi = qi + energ (srf_state(c)%tssub(i,k)-Tffresh, saltz(k))
               end do
               qi = qi/ni
               if (inout == 1) then
                  hs = snowhice(i,c) + snowfall(i,c)*dtime
               else
                  hs = snowhice(i,c)
               end if
               totenerg(i,c) = (-rLfs*hs*rhofresh/rhos + qi*sicthk(i,c))*aice(i,c)
! if (inout == 1) then
!    write(6,*)'i,c,snowterm2a=',i,c,(-rLfs*hs*rhofresh/rhos)*aice(i,c)
!    write(6,*)'i,c,iceterm2a=',i,c,qi*sicthk(i,c)*aice(i,c)
! else
!    write(6,*)'i,c,snowterm2b=',i,c,(-rLfs*hs*rhofresh/rhos)*aice(i,c)
!    write(6,*)'i,c,iceterm2b=',i,c,qi*sicthk(i,c)*aice(i,c)
! end if
            end if
         end do
         if (inout == 1) then
            do i=1,ncol
               deltae(i,c)    = totenerg(i,c)
               deltaaice(i,c) = aice(i,c)
            end do
         else
            do i=1,ncol
               deltae(i,c)    = totenerg(i,c) - deltae(i,c)
               deltaaice(i,c) = aice(i,c) - deltaaice(i,c)
               erate(i,c)     = deltae(i,c)/dtime
            end do
            call outfld ('DELTAICE', deltaaice(1,c), pcols, c)
            call outfld ('NRGICE  ', totenerg(1,c), pcols, c)
            call outfld ('IIERATE ', erate(1,c), pcols, c)
         end if
      end do

      call gather_chunk_to_field (1, 1, 1, plon, totenerg, totenerg_field)

      if (masterproc) then
         mean = 0.
         wght = 0.
         do j=1,plat
            wt = w(j)/nlon(j)
            do i=1,nlon(j)
               if (landfrac_field(i,j) < 1.) then
                  tmp = wt*(1. - landfrac_field(i,j))
                  mean = mean + tmp*totenerg_field(i,j)
                  wght = wght + tmp
               end if
            end do
         end do
         mean = mean/wght
      end if
#if ( defined SPMD )
      call mpibcast (mean, 1, mpir8, 0, mpicom)
#endif
      sea_ice_energy = mean
      call t_stopf ('sea_ice_energy')
      return
   end function sea_ice_energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine gmean_ice (srfflx, srf_state, fsns, snowfall, F_ice, &
                         F_ocn, F_frzmlt, dtime, deltae, evap, aiceinit)
!
! Arguments
!
      type(srfflx_parm  ), intent(in) :: srfflx(begchunk:endchunk)
      type(surface_state), intent(in) :: srf_state(begchunk:endchunk)
      real(r8), intent(in) :: fsns(pcols,begchunk:endchunk)
      real(r8), intent(in) :: snowfall(pcols,begchunk:endchunk)
      real(r8), intent(out) :: F_ice
      real(r8), intent(out) :: F_ocn
      real(r8), intent(out) :: F_frzmlt
      integer, intent(in) :: dtime
      real(r8), intent(in) :: deltae(pcols,begchunk:endchunk)
      real(r8), intent(in) :: evap(pcols,begchunk:endchunk)
      real(r8), intent(in) :: aiceinit(pcols,begchunk:endchunk)
!
! Local workspace
!
      real(r8) :: imbalance(pcols,begchunk:endchunk)
      real(r8) :: F_ice_chunk(pcols,begchunk:endchunk)
      real(r8) :: F_ocn_chunk(pcols,begchunk:endchunk)
      real(r8) :: F_frzmlt_chunk(pcols,begchunk:endchunk)
      real(r8) :: F_ice_field(plon,plat)
      real(r8) :: F_ocn_field(plon,plat)
      real(r8) :: F_frzmlt_field(plon,plat)
      real(r8) :: wght
      real(r8) :: wt
      real(r8) :: tmp
      real(r8) :: netflux

      integer :: i, c, j, ncol

      call t_startf ('gmean_ice')

!$OMP PARALLEL DO PRIVATE (C, I, NCOL, NETFLUX)
      do c=begchunk,endchunk
         ncol = get_ncols_p(c)
         do i=1,ncol
            F_ice_chunk(i,c) = 0.
            F_ocn_chunk(i,c) = 0.
            F_frzmlt_chunk(i,c) = 0.
            if (landfrac(i,c) < 1.) then
               F_ocn_chunk(i,c) = Focn(i,c)
               F_frzmlt_chunk(i,c) = max (0._r8, frzmlt(i,c))
            end if
            if (aiceinit(i,c) > 0.) then
!BPB note: CAM2 sign conventions on fluxes
               F_ice_chunk(i,c) = aiceinit(i,c)* &
                      (fsns(i,c) + srf_state(c)%flwds(i) - &
                       srfflx(c)%lwup(i) + evap(i,c)*Lvap - srfflx(c)%shf(i))
!               write(6,*)'i,c,aiceinit2=',i,c,aiceinit(i,c)
!               write(6,*)'i,c,aice2=',i,c,aice(i,c)
!               write(6,*)'i,c,fswabs2=',i,c,fsns(i,c)
!               write(6,*)'i,c,lwsum2=',i,c,srf_state(c)%flwds(i)-srfflx(c)%lwup(i)
!               write(6,*)'i,c,shflx2=',i,c,srfflx(c)%shf(i)
!               write(6,*)'i,c,evaplvap2=',i,c,evap(i,c)*Lvap
!               write(6,*)'i,c,focn2=',i,c,F_ocn_chunk(i,c)
!               write(6,*)'i,c,frzmlt2=',i,c,F_frzmlt_chunk(i,c)
            endif
            netflux = F_ice_chunk(i,c) - F_ocn_chunk(i,c) - F_frzmlt_chunk(i,c)
            imbalance(i,c) = deltae(i,c)/dtime - netflux
         end do
         call outfld ('F_ICE',    F_ice_chunk(1,c),    pcols, c)
         call outfld ('F_OCN',    F_ocn_chunk(1,c),    pcols, c)
         call outfld ('FRZMLTMX', F_frzmlt_chunk(1,c), pcols, c)
         call outfld ('IMBAL   ', imbalance(1,c), pcols, c)
      end do

      call gather_chunk_to_field (1, 1, 1, plon, F_ice_chunk, F_ice_field)
      call gather_chunk_to_field (1, 1, 1, plon, F_ocn_chunk, F_ocn_field)
      call gather_chunk_to_field (1, 1, 1, plon, F_frzmlt_chunk, F_frzmlt_field)

      if (masterproc) then
         F_ice    = 0.
         F_ocn    = 0.
         F_frzmlt = 0.
         wght     = 0.
         do j=1,plat
            wt = w(j)/nlon(j)
            do i=1,nlon(j)
               if (landfrac_field(i,j) < 1.) then
                  tmp = wt*(1. - landfrac_field(i,j))
                  F_ice    = F_ice    + tmp*F_ice_field(i,j)
                  F_ocn    = F_ocn    + tmp*F_ocn_field(i,j)
                  F_frzmlt = F_frzmlt + tmp*F_frzmlt_field(i,j)
                  wght     = wght     + tmp
               end if
            end do
         end do
         F_ice    = F_ice/wght
         F_ocn    = F_ocn/wght
         F_frzmlt = F_frzmlt/wght
      end if

#if ( defined SPMD )
      call mpibcast (F_ice, 1, mpir8, 0, mpicom)
      call mpibcast (F_ocn, 1, mpir8, 0, mpicom)
      call mpibcast (F_frzmlt, 1, mpir8, 0, mpicom)
#endif
      call t_stopf ('gmean_ice')
      return
   end subroutine gmean_ice
end module ice_globalcalcs
