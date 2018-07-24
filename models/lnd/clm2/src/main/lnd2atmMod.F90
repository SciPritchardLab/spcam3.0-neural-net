#include <misc.h>
#include <preproc.h>

module lnd2atmMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: tradMod
! 
! !DESCRIPTION: 
! Compute l2a component of gridcell derived type
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: lnd2atm         
!
! !REVISION HISTORY:
! 03-04-27 : Created by Mariana Vertenstein
!
!EOP
!----------------------------------------------------------------------- 

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd2atm
!
! !INTERFACE: subroutine lnd2atm(init)
  subroutine lnd2atm(init)
!
! !DESCRIPTION: 
! Compute l2a component of gridcell derived type
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use clmpoint
    use clm_varcon, only : sb
!
! !ARGUMENTS:
    implicit none
    logical, optional, intent(in) :: init  ! if true=>only set a subset of arguments
!
! !REVISION HISTORY:
! 03-04-27 : Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: pi,ci,li,gi  !indices
#if (defined DUST)
    integer :: n            !indices
#endif
    integer :: begg,endg    !indices  
    real(r8), dimension(:), pointer :: wt       ! pointer to array of pft weights
    type(gridcell_type)       , pointer :: g    ! local pointer to derived subtype
    type(landunit_type)       , pointer :: l    ! local pointer to derived subtype
    type(column_type)         , pointer :: c    ! local pointer to derived subtype
    type(lnd2atm_flux_type)   , pointer :: l2af ! local pointer to derived subtype
    type(lnd2atm_state_type)  , pointer :: l2as ! local pointer to derived subtype
    type(gridcell_pstate_type), pointer :: gps  ! local pointer to derived subtype
    type(landunit_pstate_type), pointer :: lps  ! local pointer to derived subtype
    type(column_pstate_type)  , pointer :: cps  ! local pointer to derived subtype
    type(column_wstate_type)  , pointer :: cws  ! local pointer to derived subtype
    type(pft_pstate_type)     , pointer :: pps  ! local pointer to derived subtype
    type(pft_mflux_type)      , pointer :: pmf  ! local pointer to derived subtype
    type(pft_wflux_type)      , pointer :: pwf  ! local pointer to derived subtype
    type(pft_eflux_type)      , pointer :: pef  ! local pointer to derived subtype 
    type(pft_estate_type)     , pointer :: pes  ! local pointer to derived subtype
#if (defined DUST)
    type(pft_dflux_type)      , pointer :: pdf  ! local pointer to derived subtype
#endif
! -----------------------------------------------------------------

    ! Determine per-process beginning and ending indices

    begg = grid1d%beg
    endg = grid1d%end

    ! Initialize gridcell land->atm components before grid averaging
    do gi = begg,endg
       g => gpoint(gi)%g
       l2as => g%l2as
       l2af => g%l2af
       l2as%albd(:) = 0.
       l2as%albi(:) = 0.
       l2as%t_ref2m = 0.
       l2as%q_ref2m = 0.
       l2as%h2osno = 0.
#if (defined DUST)
       l2as%fv = 0.
       l2as%ram1 = 0.
#endif
       l2af%taux = 0.
       l2af%tauy = 0. 
       l2af%eflx_lh_tot = 0.
       l2af%eflx_sh_tot = 0.
       l2af%eflx_lwrad_out = 0.
       l2af%qflx_evap_tot = 0.
       l2af%fsa = 0.
#if (defined DUST)
       l2af%flxdst(:) = 0.
#endif
       l2as%t_rad = 0.	
    end do

    ! Compute gridcell averages. Note that gridcell value for the 
    ! radiative temperature (l2as%t_rad) must be computed after the
    ! gridcell average of eflx_lwrad_out is computed.

   if (present(init)) then
      if (init) then
         write(6,*)'using intialize version of lnd2atm' 
!$OMP PARALLEL DO PRIVATE(gi,g,gps,l2as,l2af,li,l,lps,ci,c,cps,cws,pi,pps,pef)
         do gi = begg,endg
            g => gpoint(gi)%g
            gps => g%gps
            l2as => g%l2as
            l2af => g%l2af
            do li = 1, gps%nlandunits
               l => g%l(li)
               lps => l%lps
               do ci = 1, lps%ncolumns
                  c => l%c(ci)
                  cps => c%cps
                  cws => c%cws
                  l2as%h2osno = l2as%h2osno + cws%h2osno/1000. * cps%wtxy
                  do pi = 1, cps%npfts
                     pps => c%p(pi)%pps
                     pef => c%p(pi)%pef
                     l2as%albd(1) = l2as%albd(1) + pps%albd(1) * pps%wtxy  
                     l2as%albd(2) = l2as%albd(2) + pps%albd(2) * pps%wtxy 
                     l2as%albi(1) = l2as%albi(1) + pps%albi(1) * pps%wtxy 
                     l2as%albi(2) = l2as%albi(2) + pps%albi(2) * pps%wtxy 
                     l2af%eflx_lwrad_out = l2af%eflx_lwrad_out + pef%eflx_lwrad_out * pps%wtxy 
                  end do
               end do
            end do
            l2as%t_rad = sqrt(sqrt(l2af%eflx_lwrad_out/sb))
         end do
      endif
   else
#if (defined DUST)
!$OMP PARALLEL DO PRIVATE(gi,g,gps,l2as,l2af,li,l,lps,ci,c,cps,cws,pi,pps,pwf,pef,pes,pmf,pdf)
#else
!$OMP PARALLEL DO PRIVATE(gi,g,gps,l2as,l2af,li,l,lps,ci,c,cps,cws,pi,pps,pwf,pef,pes,pmf)
#endif
      do gi = begg,endg
         g => gpoint(gi)%g
         gps => g%gps
         l2as => g%l2as
         l2af => g%l2af
         do li = 1, gps%nlandunits
            l => g%l(li)
            lps => l%lps
            do ci = 1, lps%ncolumns
               c => l%c(ci)
               cps => c%cps
               cws => c%cws
               l2as%h2osno = l2as%h2osno + cws%h2osno/1000. * cps%wtxy
                do pi = 1, cps%npfts
                  pps => c%p(pi)%pps
                  pwf => c%p(pi)%pwf
                  pef => c%p(pi)%pef
                  pes => c%p(pi)%pes
                  pmf => c%p(pi)%pmf
#if (defined DUST)
                  pdf => c%p(pi)%pdf
#endif
                  l2as%albd(1) = l2as%albd(1) + pps%albd(1) * pps%wtxy  
                  l2as%albd(2) = l2as%albd(2) + pps%albd(2) * pps%wtxy 
                  l2as%albi(1) = l2as%albi(1) + pps%albi(1) * pps%wtxy 
                  l2as%albi(2) = l2as%albi(2) + pps%albi(2) * pps%wtxy 
                  l2as%t_ref2m = l2as%t_ref2m + pes%t_ref2m * pps%wtxy 
                  l2as%q_ref2m = l2as%q_ref2m + pes%q_ref2m * pps%wtxy 
                  l2af%taux = l2af%taux + pmf%taux * pps%wtxy 
                  l2af%tauy = l2af%tauy + pmf%tauy * pps%wtxy 
                  l2af%eflx_lh_tot = l2af%eflx_lh_tot + pef%eflx_lh_tot * pps%wtxy 
                  l2af%eflx_sh_tot = l2af%eflx_sh_tot + pef%eflx_sh_tot * pps%wtxy 
                  l2af%eflx_lwrad_out = l2af%eflx_lwrad_out + pef%eflx_lwrad_out * pps%wtxy 
                  l2af%qflx_evap_tot  = l2af%qflx_evap_tot  + pwf%qflx_evap_tot  * pps%wtxy 
                  l2af%fsa = l2af%fsa + pef%fsa * pps%wtxy 
#if (defined DUST)
! add for dust (nmm)
                  l2as%fv=l2as%fv+ pps%fv * pps%wtxy 
                  l2as%ram1=l2as%ram1+ pps%ram1 * pps%wtxy 
                  do n=1,ndst
                     l2af%flxdst(n)=l2af%flxdst(n)+ pdf%flx_mss_vrt_dst(n) * pps%wtxy 
                  enddo
#endif
               end do
            end do
         end do
         l2as%t_rad = sqrt(sqrt(l2af%eflx_lwrad_out/sb))
      end do
    endif

  end subroutine lnd2atm

end module lnd2atmMod
