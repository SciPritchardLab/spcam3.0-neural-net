c CVS: $Id: ice_diagnostics.F,v 1.10.6.1 2003/02/27 00:58:49 rosinski Exp $
c CVS: $Source: /fs/cgd/csm/models/CVS.REPOS/ice/csim4_src/source/ice_diagnostics.F,v $
c CVS: $Name:  $
!
!BOP
!
! !MODULE: ice_diagnostics - diagnostic information output during run
!
! !DESCRIPTION:
!
! Diagnostic information output during run
!
! !REVISION HISTORY:
!
! authors: Elizabeth C. Hunke, LANL 
!          Bruce P. Briegleb, NCAR
!
! !INTERFACE:
!
      module ice_diagnostics
!
! !USES:
!
      use ice_constants
!
!EOP
!
      implicit none

      ! for water and heat budgets
      real (kind=dbl_kind) ::
     &   totmn   ! total ice/snow water mass (nh)
     &,  totms   ! total ice/snow water mass (sh)
     &,  toten   ! total ice/snow enthalpy (J)
     &,  totes   ! total ice/snow enthalpy (J)

      ! printing info for routine print_state
      character (len=20) :: plabel
      integer, parameter ::
     &  check_step = 100000
     &, ip = 51
     &, jp = 37
     &, mtask = 2
!
      contains
!
!BOP
!
! !IROUTINE: print_state - print ice state for specified grid point
!
! !INTERFACE:
!
      subroutine print_state(plabel,i,c)
!
! !DESCRIPTION:
!
! This routine is useful for debugging call to it should be inserted 
! in the form (after thermo, for example)
!
!      do j=jlo,jhi    \\
!      -do i=ilo,ihi   \\ 
!      ---call to\_column(i,j)    \\
!      ---plabel = 'post thermo' \\
!      ---if (istep1.ge.check\_step.and.i.eq.ip.and.j.eq.jp  \\
!      ------.and.my\_task.eq.mtask)  \\
!      ---call print\_state(plabel,i,j) \\
!      -enddo  \\
!      enddo   \\
!
! 'use ice\_diagnostics' may need to be inserted also, and \\
! 'use ice\_calendar' if it is not already being used
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_kinds_mod
      use comsrf, only: aice, sicthk, snowhice, surface_state2d, 
     &                  frzmlt, tsocn
      use phys_grid, only: get_lat_p, get_lon_p
      use time_manager, only: get_nstep
!
! !INPUT/OUTPUT PARAMETERS:
!
      character (len=20), intent(in) :: plabel
      integer, intent(in) :: i,c
!
!EOP
!
      real (kind=dbl_kind) :: edebug

      integer :: nc, layer
      integer :: nunit
      integer :: lon, lat
      real (kind=dbl_kind) :: z

      lat = get_lat_p (c,i)
      lon = get_lon_p (c,i)

      nunit = 6

      write (nunit,*) plabel,' nstep = ', get_nstep ()
      write (nunit,*) 'lon,lat=',lon,lat
      write (nunit,*) 'aice ',aice(i,c)
!      write (nunit,*) 'tmask =',tmask(i,j)

      write (nunit,*) 'nc=',1
      write (nunit,*) 'ain ',aice(i,c)
      write (nunit,*) 'vin ',aice(i,c)*sicthk(i,c)
      write (nunit,*) 'hin ',sicthk(i,c)
      write (nunit,*) 'vsn ',aice(i,c)*snowhice(i,c)
      write (nunit,*) 'hsn ',snowhice(i,c)
!      write (nunit,*) 'Tsfn',Tsfn(nc)
      write (nunit,*) ' '

!      edebug = c0
!      do layer = 1,ni
!        write (nunit,*) 'ein, cat ',nc,' layer ',layer,ein(layer,nc)
!        edebug = edebug+ein(layer,nc)
!      enddo
!      write (nunit,*) 'eice(i,j)',edebug 

      write (nunit,*) ' '
      write (nunit,*) 'atm states and fluxes'
      write (nunit,*) '            zlvl  = ',surface_state2d(c)%zbot(i)
!      write (nunit,*) '            uatm  = ',uatm(i,c)
!      write (nunit,*) '            vatm  = ',vatm(i,c)
      write (nunit,*) '            potT  = ',surface_state2d(c)%thbot(i)
      write (nunit,*) '            Tair  = ',surface_state2d(c)%tbot(i)
      write (nunit,*) '            Qa    = ',surface_state2d(c)%qbot(i)
!      write (nunit,*) '            rhoa  = ',rhoa(i,c)
      write (nunit,*) '            swvdr = ',surface_state2d(c)%sols(i)
      write (nunit,*) '            swvdf = ',surface_state2d(c)%solsd(i)
      write (nunit,*) '            swidr = ',surface_state2d(c)%soll(i)
      write (nunit,*) '            swidf = ',surface_state2d(c)%solld(i)
      write (nunit,*) '            Flw   = ',surface_state2d(c)%flwds(i)
!      write (nunit,*) '            Frain = ',Frain(i,c)
      z = surface_state2d(c)%precsc(i) + surface_state2d(c)%precsl(i)
      write (nunit,*) '            Fsnow = ',z
      write (nunit,*) ' '
      write (nunit,*) 'ocn states and fluxes'
      write (nunit,*) '            frzmlt = ',frzmlt(i,c)
      write (nunit,*) '            sst    = ',tsocn(i,c)
!      write (nunit,*) '            sss    = ',sss(i,c)
      write (nunit,*) '            Tf     = ',-1.8
      write (nunit,*) ' '
      write (nunit,*) 'srf states and fluxes'
!      write (nunit,*) '            Tref      = ',Tref(i,c)
!      write (nunit,*) '            Fsensible = ',Fsensible(i,c)
!      write (nunit,*) '            Flatent   = ',Flatent(i,c)
!      write (nunit,*) '            evap      = ',evap(i,c)
!      write (nunit,*) '            Flwout    = ',Flwout(i,c)
      write (nunit,*) ' '

      end subroutine print_state
!
      end module ice_diagnostics
!
