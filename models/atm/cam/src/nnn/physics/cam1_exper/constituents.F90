#include <misc.h>
#include <params.h>

module constituents

!----------------------------------------------------------------------- 
! 
! Purpose: Contains data and functions for manipulating advected and non-advected constituents
!
! Author: B.A. Boville
! 
!-----------------------------------------------------------------------
  use shr_kind_mod, only: r8 => shr_kind_r8
  use physconst, only: r_universal
  use pmgrid,    only: masterproc

  implicit none

  private
  save
!
! Public interfaces
!
  public cnst_add          ! add a constituent to the list of advected (or nonadvected) constituents
  public cnst_get_ind      ! get the index of a constituent
  public cnst_read_iv      ! query whether constituent initial values are read from initial file
  public cnst_chk_dim      ! check that number of constituents added equals dimensions (pcnst, pnats)

  integer, parameter, public :: pcnst  = PCNST      ! number of advected constituents (including water vapor)
  integer, parameter, public :: pnats  = PNATS      ! number of non-advected constituents
  integer, parameter, public :: ppcnst = pcnst+pnats! total number of constituents

  integer, parameter, public :: advected = 0        ! type value for constituents which are advected
  integer, parameter, public :: nonadvec = 1        ! type value for constituents which are not advected

  character(len=8),public :: cnst_name(ppcnst)      ! constituent names (including any non-advected)
  character(len=128),public :: cnst_longname(ppcnst)! long name of constituents

! Public data

! Namelist variables
  logical, public :: readtrace = .true.             ! true => obtain initial tracer data from IC file

!
! Constants for each tracer
  real(r8), public :: cnst_cp  (ppcnst)             ! specific heat at constant pressure (J/kg/K)
  real(r8), public :: cnst_cv  (ppcnst)             ! specific heat at constant volume (J/kg/K)
  real(r8), public :: cnst_mw  (ppcnst)             ! molecular weight (kg/kmole)
  real(r8), public :: cnst_rgas(ppcnst)             ! gas constant ()
  real(r8), public :: qmin     (ppcnst)             ! minimum permitted constituent concentration (kg/kg)
  real(r8), public :: qmincg   (ppcnst)             ! for backward compatibility only

!++bee - temporary...
! Lists of tracer names and diagnostics
   character(len=8), public :: hadvnam(pcnst)        ! names of horizontal advection tendencies
   character(len=8), public :: vadvnam(pcnst)        ! names of vertical advection tendencies
   character(len=8), public :: dcconnam(ppcnst)      ! names of convection tendencies
   character(len=8), public :: fixcnam(pcnst)        ! names of species slt fixer tendencies
   character(len=8), public :: tendnam(pcnst)        ! names of total tendencies of species
   character(len=8), public :: sflxnam(ppcnst)       ! names of surface fluxes of species
   character(len=8), public :: tottnam(pcnst)        ! names for horz + vert + fixer tendencies
   character(len=8), public :: qphystendnam(ppcnst)   ! names for physics tendencies
!--bee

! Private data

  integer :: padv = 0                      ! index pointer to last advected tracer
  integer :: pnad = pcnst                  ! index pointer to last non-advected tracer
  logical :: read_init_vals(ppcnst)        ! true => read initial values from initial file

CONTAINS

!==============================================================================
  subroutine cnst_add (name, type, mwc, cpc, qminc, ind, longname, readiv)
!----------------------------------------------------------------------- 
! 
! Purpose: Register a constituent for treament by physical parameterizations and 
!          transport (if type=advected)
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author:  B.A. Boville
! 
!-----------------------------Arguments---------------------------------
!
    character(len=*), intent(in) :: &
       name      ! constituent name used as variable name in history file output (8 char max)
    character(len=*), intent(in), optional :: &
       longname  ! value for long_name attribute in netcdf output (128 char max, defaults to name)
    logical,          intent(in), optional :: &
       readiv    ! true => read initial values from initial file (default: true)

    integer, intent(in)    :: type   ! flag indicating advected or nonadvected
    real(r8),intent(in)    :: mwc    ! constituent molecular weight (kg/kmol)
    real(r8),intent(in)    :: cpc    ! constituent specific heat at constant pressure (J/kg/K)
    real(r8),intent(in)    :: qminc  ! minimum value of mass mixing ratio (kg/kg)
!                                        normally 0., except water 1.E-12, for radiation.

    integer, intent(out)   :: ind    ! global constituent index (in q array)

!-----------------------------------------------------------------------

! set tracer index and check validity, advected tracer
    if (type == advected) then
       padv = padv+1
       ind  = padv
       if (padv > pcnst) then
          write(6,*) 'CNST_ADD: advected tracer index greater than pcnst = ', pcnst
          call endrun
       end if

! set tracer index and check validity, non-advected tracer
    else if (type == nonadvec) then
       pnad = pnad+1
       ind  = pnad
       if (pnad > ppcnst) then
          write(6,*) 'CNST_ADD: non-advected tracer index greater than pcnst+pnats = ', ppcnst
          call endrun
       end if

! unrecognized type value
    else
       write(6,*) 'CNST_ADD, input type flag not valid, type=', type
       call endrun
    end if

! set tracer name and constants
    cnst_name(ind) = name
    if ( present(longname) )then
       cnst_longname(ind) = longname
    else
       cnst_longname(ind) = name
    end if

    read_init_vals(ind) = readtrace
    if ( present(readiv) ) read_init_vals(ind) = readiv

    cnst_cp  (ind) = cpc
    cnst_mw  (ind) = mwc
    qmin     (ind) = qminc
    qmincg   (ind) = qminc
    if (ind == 1) qmincg = 0.  ! This crap is replicate what was there before ****

    cnst_rgas(ind) = r_universal * mwc
    cnst_cv  (ind) = cpc - cnst_rgas(ind)

    return
  end subroutine cnst_add

!==============================================================================
  subroutine cnst_get_ind (name, ind)
!----------------------------------------------------------------------- 
! 
! Purpose: Get the index of a constituent 
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author:  B.A. Boville
! 
!-----------------------------Arguments---------------------------------
!
    character(len=*), intent(in) :: name ! constituent name

    integer, intent(out)   :: ind    ! global constituent index (in q array)

!---------------------------Local workspace-----------------------------
    integer :: m                                   ! tracer index

!-----------------------------------------------------------------------

! Find tracer name in list
    do m = 1, ppcnst
       if (name == cnst_name(m)) then
          ind  = m
          return
       end if
    end do

! Unrecognized name
    write(6,*) 'CNST_GET_IND, name:', name,  ' not found in list:', cnst_name(:)
    call endrun

  end subroutine cnst_get_ind

!==============================================================================
  function cnst_read_iv(m)
!----------------------------------------------------------------------- 
! 
! Purpose: Query whether constituent initial values are read from initial file.
! 
! Author:  B. Eaton
! 
!-----------------------------Arguments---------------------------------
!
    integer, intent(in) :: m    ! constituent index

    logical :: cnst_read_iv     ! true => read initial values from inital file
!-----------------------------------------------------------------------

    cnst_read_iv = read_init_vals(m)
 end function cnst_read_iv

!==============================================================================
  subroutine cnst_chk_dim
!----------------------------------------------------------------------- 
! 
! Purpose: Check that the number of registered constituents of each type is the
!          same as the dimension
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author:  B.A. Boville
! 
    integer i
!-----------------------------------------------------------------------
!
    if (padv /= pcnst) then
       write(6,*)'CNST_CHK_DIM: number of advected tracer ',padv, ' not equal to pcnst = ',pcnst
       call endrun ()
    endif
    if (pnad /= ppcnst) then
       write(6,*)'CNST_CHK_DIM: number of non-advected tracers ',pnad, ' not equal to pcnst+pnats = ', &
            ppcnst
       call endrun ()
    endif

    if (masterproc) then
       print *, 'Advected constituent list:'
       do i = 1, pcnst
          write(6,'(i4,2x,a8,2x,a128)') i, cnst_name(i), cnst_longname(i)
       end do
       if (ppcnst > pcnst) then
          print *, 'Nonadvected constituent list:'
          do i = pcnst+1, ppcnst
             write(6,'(i4,2x,a8,2x,a128)') i, cnst_name(i), cnst_longname(i)
          end do
       else
          print *, 'No nonadvected constituents'
       end if
    end if
  end subroutine cnst_chk_dim

end module constituents
