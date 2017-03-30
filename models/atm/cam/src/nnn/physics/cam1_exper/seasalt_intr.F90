#include <misc.h>
#include <params.h>

module seasalt_intr

!---------------------------------------------------------------------------------
! Module to interface the aerosol parameterizations with CAM
! written by PJR (extensively modified from chemistry module)
!---------------------------------------------------------------------------------

  use shr_kind_mod,only: r8 => shr_kind_r8
  use pmgrid,      only: plat, plev, plevp, masterproc
  use ppgrid,      only: pcols, pver
  use physconst,   only: mwdry, mwh2o
  use constituents,only: ppcnst, cnst_add, cnst_name, advected, cnst_get_ind
    
  implicit none

  private          ! Make default type private to the module

  save

! Public interfaces
!
  public seasalt_register_cnst                        ! register consituents
  public seasalt_implements_cnst                      ! returns true if consituent is implemented by this package
  public seasalt_init_cnst                            ! initialize mixing ratios if not read from initial file
  public seasalt_initialize                           ! initialize (history) variables
  public seasalt_srcsnk                               ! production and loss of seasalt

  integer ixsslt
  integer, parameter :: ncnst=1                      ! number of constituents
  character(len=8), dimension(ncnst), parameter :: & ! constituent names
     cnst_names = (/'SSLT'/)

  real(r8)  gravit

contains

!===============================================================================
  subroutine seasalt_register_cnst
!----------------------------------------------------------------------- 
! 
! Purpose: register seasalt
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: P. J. Rasch
! 
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
    integer :: m                                   ! tracer index
    real(r8), parameter :: one  = 1._r8
    real(r8), parameter :: zero  = 0._r8
!-----------------------------------------------------------------------

! Set names of variables undergoing evolution
    call cnst_add(cnst_names(1), advected, one, one, zero, ixsslt) 
    write (6,*) ' seasalt_register_cnst: start index for seasalt is ', ixsslt

    return
  end subroutine seasalt_register_cnst


!=======================================================================
  function seasalt_implements_cnst(name)
!----------------------------------------------------------------------- 
! 
! Purpose: return true if specified constituent is implemented by this 
!          package
! 
! Author: T. Henderson
! 
!-----------------------------------------------------------------------
     implicit none
!-----------------------------Arguments---------------------------------

     character(len=*), intent(in) :: name   ! constituent name
     logical :: seasalt_implements_cnst     ! return value
!---------------------------Local workspace-----------------------------
     integer :: m
!-----------------------------------------------------------------------

     seasalt_implements_cnst = .false.
     do m = 1, ncnst
        if (name == cnst_names(m)) then
           seasalt_implements_cnst = .true.
           return
        end if
     end do
  end function seasalt_implements_cnst


!=======================================================================
  subroutine seasalt_init_cnst(name, q)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Set initial mass mixing ratios.
!
!-----------------------------------------------------------------------
    implicit none
!-----------------------------Arguments---------------------------------
    
    character(len=*), intent(in) :: name         ! constituent name
    
    real(r8), intent(out) :: q(:,:,:)            !  mass mixing ratio
!-----------------------------------------------------------------------
    
    if ( name == cnst_names(1) ) then
       q = 0._r8
    end if

  end subroutine seasalt_init_cnst


!===============================================================================
  subroutine seasalt_initialize 
!----------------------------------------------------------------------- 
! 
! Purpose: initialize parameterization of seasalt 
!          (declare history variables)
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: NCAR CMS
! 
!-----------------------------------------------------------------------
    use history,    only: addfld, add_default, phys_decomp
    use surface,    only: inisflx
    use shr_const_mod,    only: SHR_CONST_RDAIR, SHR_CONST_G, SHR_CONST_REARTH
    use time_manager, only: get_curr_date, get_perp_date, get_curr_calday, &
                            is_perpetual
    use commap, only: clat, clon, w
    use ppgrid, only: pcols, begchunk, endchunk


    implicit none

!---------------------------Local workspace-----------------------------
    integer :: m                               ! tracer index

    gravit = SHR_CONST_G
    write (6,*) ' seasalt_initialize: ', pcols

  end subroutine seasalt_initialize


!===============================================================================
  subroutine seasalt_srcsnk (state, dt, u10, ocnfrac, ptend)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Interface to emission of all seasalts
! 
! Method: 
! 
! Author: Phil Rasch
! 
!-----------------------------------------------------------------------
    use physics_types, only: physics_state, physics_ptend
    use phys_grid,     only: get_lon_all_p, get_lat_all_p, get_rlat_all_p
!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
!
! Arguments:
!
    type(physics_state), intent(in ) :: state          ! Physics state variables
    real(r8),            intent(in)  :: dt             ! time step
    real(r8),            intent(in)  :: u10(pcols)     ! 10 meter wind (m/s)
    real(r8),            intent(in)  :: ocnfrac(pcols) ! ocean fraction
    type(physics_ptend), intent(inout) :: ptend        ! indivdual parameterization tendencies

! local
    integer i
    integer k
    integer kstar
    integer lchnk
    integer ncol
    real(r8) atmp(pcols)
    real(r8) alpha
    real(r8) beta
    real(r8) dz
    real(r8) gamma
    real(r8) h
    real(r8) sslt(pcols,pver) ! prescribed sea salt profile
    real(r8) u
    real(r8) z
    real(r8) zib
    real(r8) zmin
!
    lchnk = state%lchnk
    ncol = state%ncol
    
    sslt = 0.

!     note that our burden is only from the sfc to about 500m 
!          (for the NCEP grid)
!     many explicit models find significant salt burdens above this level
!     to maybe 5km, although the studies of blanchard & woodcock show
!     concentrations at least an order of magnitude smaller than the
!     500m level values, and their fit is identified as being 
!          appropriate for 1 to 300m.
!

      do k = pver,1,-1
         zmin = 1.e36
!
!        write (6,*) ' points with salt, ug/m3, zt, zb, ug/m2 at lev ', k
!
         alpha = 5.e-6          ! gm/m3
         beta = 6.3e-6          ! /m
         do i = 1,ncol
            z = max(min(300._r8,state%zm(i,k)),1._r8)
            zmin = min(zmin,state%zm(i,k))
            u = max(min(14._r8,u10(i)),1._r8)
            dz = state%zi(i,k)-state%zi(i,k+1)
            gamma = 0.21-0.39*log10(u)
!
!           the next line follows the profile of blanchard and woodcock
!           sslt(i,k) = alpha*(beta*z)**gamma ! grams/m3
!           write (6,*) ' sslt midpoint (micrograms/m3)', 
!    $                  sslt(i,k)*1.e6, z
!            
            zib = max(state%zi(i,k+1),1._r8)
!
!           the next line is the mean value over the model layers 
!           of the profile of blanchard and woodcock
!
            sslt(i,k) = alpha/(dz*(gamma+1.)*beta) &
                 *( (beta*state%zi(i,k))**(gamma+1.) &
                   -(beta*zib)**(gamma+1.) &
                  )*ocnfrac(i)
         end do
!
!        write (6,*) ' zmin for level k ', zmin, k
!   
         if (zmin.gt.300) go to 100  ! if all layers are 300 meters at this eta level, break
      end do
 100  continue
      
!
!        above this level we just make the sea salt get small fast up to 5 km
!        then forget it
!
         h = 500.      ! 500 meter efolding depth (sort of from fig 7 from b&w)
         kstar = k
         do i = 1,ncol
            atmp(i) = sslt(i,kstar)*exp(state%zm(i,kstar)/h)
         end do
         do k = kstar-1,1,-1
            zmin = 1.e36
            do i = 1,ncol
               zmin = min(zmin,state%zm(i,k))
               sslt(i,k) = atmp(i)*exp(-state%zm(i,k)/h)*ocnfrac(i)
            end do
            if (zmin.gt.5.e3) go to 200
         end do
 200     continue
      
      
!++bee, convert sslt concentration from g/m3 --> kg/kg
      do k = 1,pver
         do i = 1,ncol
            dz = state%zi(i,k)-state%zi(i,k+1)
            sslt(i,k) = sslt(i,k)*.001*dz*state%rpdel(i,k)*gravit
         end do
      end do
!--bee

      ptend%q(:ncol,:,ixsslt) = (sslt(:ncol,:)-state%q(:ncol,:,ixsslt))/dt

  end subroutine seasalt_srcsnk 

end module seasalt_intr
