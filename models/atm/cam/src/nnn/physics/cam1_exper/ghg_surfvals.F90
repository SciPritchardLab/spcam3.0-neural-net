module ghg_surfvals

!-----------------------------------------------------------------------------------
! Purpose: Provides greenhouse gas (ghg) values at the Earth's surface.
!          These values may be time dependent.
!
! Author: Brian Eaton (assembled module from existing scattered code pieces)
!-----------------------------------------------------------------------------------

   use shr_kind_mod, only: r8=>shr_kind_r8
   use pmgrid,       only: masterproc
   use physconst,    only: mwdry, mwco2

   implicit none
   private
   save

! Public methods

   public ::&
      ghg_surfvals_init,  &! initialize options that depend on namelist input
      ghg_surfvals_set,   &! set ghg surface values when ramp option is on
      ghg_surfvals_ramp    ! returns true when ramp option is on

! Public data

   real(r8), public :: co2vmr = 3.550e-4         ! co2   volume mixing ratio 
   real(r8), public :: co2mmr                    ! co2   mass mixing ratio 
   real(r8), public :: n2ovmr = 0.311e-6         ! n2o   volume mixing ratio 
   real(r8), public :: ch4vmr = 1.714e-6         ! ch4   volume mixing ratio 
   real(r8), public :: f11vmr = 0.280e-9         ! cfc11 volume mixing ratio 
   real(r8), public :: f12vmr = 0.503e-9         ! cfc12 volume mixing ratio 

   character(len=16), public :: scenario_ghg = 'FIXED' ! 'FIXED' or 'RAMPED'
   integer, public  :: rampYear_ghg = 0          ! ramped gases fixed at this year (if > 0)

! Private module data

   logical :: doRamp_ghg    ! true => turn on ramping for ghg
   logical :: fixYear_ghg   ! true => Ramped gases fixed at specified year.
   real(r8), parameter :: rmwco2 = mwco2/mwdry   ! ratio of molecular weights of co2 to dry air

#include <ramp_ghg_bau.h>

!=========================================================================================
contains
!=========================================================================================

subroutine ghg_surfvals_init()

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize the ramp options that are controlled by namelist input.
! Set surface values at initial time.
! N.B. This routine must be called after the time manager has been initialized
!      since ghg_surfvals_set calls time manager methods.
! 
! Author: B. Eaton - merged code from parse_namelist and rampnl_ghg.
! 
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------

   if (scenario_ghg == 'FIXED') then
      doRamp_ghg = .false.
      if (masterproc) &
         write(6,*)'ghg_surfvals_init: ghg surface values are fixed as follows'
   else if (scenario_ghg == 'RAMPED') then
      doRamp_ghg = .true.
      if (masterproc) &
         write(6,*)'ghg_surfvals_init: ghg surface values are interpolated from scenario: ',scenario
   else
      if (masterproc) &
         write(6,*)'ghg_surfvals_init: input namelist SCENARIO_GHG must be set to either FIXED or RAMPED'
      call endrun
   endif

   if (doRamp_ghg) then
      fixYear_ghg = .false.
      if ( rampYear_ghg > 0 ) then
         fixYear_ghg = .true.
         if (masterproc) &
            write(6,*) '  FIXED values from year ',rampYear_ghg
      else
         if (masterproc) &
            write(6,*) '  RAMPED values initialized to'
      end if
      ! Initialize surface values when using ramp option.
      call ghg_surfvals_set()
   else
      co2mmr = rmwco2 * co2vmr  ! co2   mass mixing ratio 
   endif

   if (masterproc) then
      write(6,*) '  co2 volume mixing ratio = ',co2vmr
      write(6,*) '  ch4 volume mixing ratio = ',ch4vmr
      write(6,*) '  n2o volume mixing ratio = ',n2ovmr
      write(6,*) '  f11 volume mixing ratio = ',f11vmr
      write(6,*) '  f12 volume mixing ratio = ',f12vmr
   end if

end subroutine ghg_surfvals_init

!=========================================================================================

function ghg_surfvals_ramp()
   implicit none
   logical :: ghg_surfvals_ramp
   ghg_surfvals_ramp = doRamp_ghg
end function ghg_surfvals_ramp

!=========================================================================================

subroutine ghg_surfvals_set()
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Computes greenhouse gas volume mixing ratios via interpolation of
! yearly input data.
! 
! Author: B. Eaton - updated ramp_ghg for use in ghg_surfvals module
! 
!-----------------------------------------------------------------------
   use time_manager, only: get_curr_date, get_curr_calday

   implicit none

!---------------------------Local variables-----------------------------

   integer yrmodel           ! model year
   integer nyrm              ! year index
   integer nyrp              ! year index
   integer :: yr, mon, day   ! components of a date
   integer :: ncdate         ! current date in integer format [yyyymmdd]
   integer :: ncsec          ! current time of day [seconds]

   real(r8) :: calday            ! current calendar day
   real(r8) doymodel             ! model day of year
   real(r8) doydatam             ! day of year for input data yrdata(nyrm)
   real(r8) doydatap             ! day or year for input data yrdata(nyrp)
   real(r8) deltat               ! delta time
   real(r8) fact1, fact2         ! time interpolation factors
   real(r8) cfcscl               ! cfc scale factor for f11

!
! ---------------------------------------------------------------------
!
   calday = get_curr_calday()
   call get_curr_date(yr, mon, day, ncsec)
   ncdate = yr*10000 + mon*100 + day
!
! determine index into input data
!
   if ( fixYear_ghg ) then
      yrmodel  = rampYear_ghg
   else
      yrmodel  = ncdate/10000
   end if

   nyrm       = yrmodel - yrdata(1) + 1
   nyrp       = nyrm + 1
!
! if current date is before yrdata(1), quit
!
   if (nyrm < 1) then
      if ( masterproc )then
         write(6,*)'ghg_surfvals_set: data time index is out of bounds'
         write(6,*)'nyrm = ',nyrm,' nyrp= ',nyrp, ' ncdate= ', ncdate
      end if
      call endrun
   endif
!
! if current date later than yrdata(ntim), call endrun.
! if want to use ntim values - uncomment the following lines
! below and comment the call to endrun and previous write
!
   if (nyrp > ntim) then
      if ( masterproc )then
         write(6,*)'ghg_surfvals_set: error - current date is past the end of ', &
                ' valid data'
      end if
      call endrun
!         write(6,*)'ghg_surfvals_set: using ghg data for ',yrdata(ntim)
!         co2vmr = co2(ntim)*1.e-06
!         ch4vmr = ch4(ntim)*1.e-09
!         n2ovmr = n2o(ntim)*1.e-09
!         f11vmr = f11(ntim)*1.e-12*(1.+cfcscl)
!         f12vmr = f12(ntim)*1.e-12
!         co2mmr = rmwco2 * co2vmr
!         return
   endif
!
! determine time interpolation factors, check sanity
! of interpolation factors to within 32-bit roundoff
! assume that day of year is 1 for all input data
!
   doymodel = yrmodel*365.    + calday
   doydatam = yrdata(nyrm)*365. + 1.
   doydatap = yrdata(nyrp)*365. + 1.
   deltat   = doydatap - doydatam
   fact1    = (doydatap - doymodel)/deltat
   fact2    = (doymodel - doydatam)/deltat

   if (abs(fact1+fact2-1.) > 1.e-6 .or. &
       fact1 > 1.000001 .or. &
       fact1 < -1.e-6 .or. &
       fact2 > 1.000001 .or. &
       fact2 < -1.e-6) then
      if ( masterproc ) &
         write(6,*)'ghg_surfvals_set: Bad fact1 and/or fact2=',fact1,fact2
      call endrun
   end if
!
! do time interpolation:
!   co2     in ppmv
!   n2o,ch4 in ppbv
!   f11,f12 in pptv
!
   co2vmr = (co2(nyrm)*fact1 + co2(nyrp)*fact2)*1.e-06
   ch4vmr = (ch4(nyrm)*fact1 + ch4(nyrp)*fact2)*1.e-09
   n2ovmr = (n2o(nyrm)*fact1 + n2o(nyrp)*fact2)*1.e-09

   cfcscl = (adj(nyrm)*fact1 + adj(nyrp)*fact2)
   f11vmr = (f11(nyrm)*fact1 + f11(nyrp)*fact2)*1.e-12*(1.+cfcscl)
   f12vmr = (f12(nyrm)*fact1 + f12(nyrp)*fact2)*1.e-12

   co2mmr = rmwco2 * co2vmr
!
! output statements of ramping
!
   if ( masterproc )then
      write(6,'(a,f8.2,6(1pe22.14))') 'calday1 = ',calday,co2vmr/1.e-06,ch4vmr/1.e-09, &
                                      n2ovmr/1.e-09
      write(6,'(a,f8.2,6(1pe22.14))') 'calday2 = ',calday,cfcscl, &
                                      (f11(nyrm)*fact1 + f11(nyrp)*fact2),f12vmr/1.e-12
   end if
   return
end subroutine ghg_surfvals_set

!=========================================================================================

end module ghg_surfvals
