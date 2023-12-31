#include <misc.h>
#include <params.h>
!----------------------------------------------------------------------- 
!
! BOP
!
! !MODULE: ice_data
!
! !DESCRIPTION:	Module to handle dealing with the ICE datasets.
!
! Public interfaces:
!
!	iceini -- Initialization and reading of dataset.
!	iceint -- Interpolate dataset ICE to current time.
!	icean --- Apply the interpolated ICE to the model state.
!
!----------------------------------------------------------------------- 

module ice_data
!
! USES:
!
  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid,    only: plon, plat, masterproc
  use ppgrid,    only: pcols, begchunk, endchunk
  use phys_grid, only: scatter_field_to_chunk, get_ncols_p
  use comsrf,    only: plevmx, icefrac,previcefrac,update_ocnice,aice
  use physconst, only: tmelt
  use commap,    only: clat, clon

  implicit none
!----------------------------------------------------------------------- 
! PUBLIC: Make default data and interfaces private
!----------------------------------------------------------------------- 
!
! ! PUBLIC MEMBER FUNCTIONS:
!
!  public icean    ! Set the surface temperature, oro, and sea-ice fraction
  public iceini   ! Initialization
  public iceint   ! Time interpolation of ICE data
!===============================================================================
!EOP
!===============================================================================
!----------------------------------------------------------------------- 
! PRIVATE: Everthing else is private to this module
!----------------------------------------------------------------------- 
  private   ! By default all data is private to this module
  integer, parameter :: toticesz=9612
  real(r8), parameter :: daysperyear = 365.0  ! Number of days in a year

  real(r8), allocatable, dimension(:,:,:) :: &
      icebdy         ! ICE values on boundary dataset (pcols,begchunk:endchunk,2)
  real(r8), allocatable, dimension(:,:) :: &
      ice            ! Interpolated model ice values (pcols,begchunk:endchunk)

  real(r8) :: cdayicem   ! Calendar day for prv. month ICE values read in
  real(r8) :: cdayicep   ! Calendar day for nxt. month ICE values read in
      
  integer :: nm,np   ! Array indices for prv., nxt month ice data
  integer :: iceid   ! netcdf id for ice variable
  integer :: lonsiz  ! size of longitude dimension on ice dataset
  integer :: levsiz  ! size of level dimension on ice dataset
  integer :: latsiz  ! size of latitude dimension on ice dataset
  integer :: timesiz ! size of time dimension on ice dataset
  integer :: np1     ! current forward time index of ice dataset
  integer :: date_ice(toticesz)! Date on ice dataset (YYYYMMDD)
  integer :: sec_ice(toticesz) ! seconds of date on ice dataset (0-86399)

  real(r8), parameter :: tsice = -1.7999 ! Freezing point of sea ice degrees C
                                         ! Use this with global ice data
 
!===============================================================================
CONTAINS
!===============================================================================

!======================================================================
! PUBLIC ROUTINES: Following routines are publically accessable
!======================================================================
!----------------------------------------------------------------------- 
! 
! BOP
!
! !IROUTINE: iceini
!
! !DESCRIPTION:
!
! Initialize the procedure for specifying sea surface temperatures
! Do initial read of time-varying ice boundary dataset, reading two
! consecutive months on either side of the current model date.
!
! Method: 
! 
! Author: L.Bath
! 
!-----------------------------------------------------------------------
!
! !INTERFACE
!
subroutine iceini
!
! !USES:
!
  use rgrid, only: nlon
  use error_messages, only: alloc_err, handle_ncerr
  use time_manager, only: get_curr_date, get_curr_calday, get_step_size, &
                          is_first_step

#if ( defined SPMD )
  use mpishorthand, only: mpicom, mpiint, mpir8
#endif
!
! EOP
!
!---------------------------Common blocks-------------------------------
#include <comctl.h>
#include <comlun.h>
!---------------------------Local variables-----------------------------
  integer dtime                 ! timestep size [seconds]
  integer dateid                ! netcdf id for date variable
  integer secid                 ! netcdf id for seconds variable
  integer londimid              ! netcdf id for longitude variable
  integer latdimid              ! netcdf id for latitude variable
  integer lonid                 ! netcdf id for longitude variable
  integer latid                 ! netcdf id for latitude variable
  integer timeid                ! netcdf id for time variable
  integer nlonid                ! netcdf id for nlon variable (rgrid)
  integer cnt3(3)               ! array of counts for each dimension
  integer strt3(3)              ! array of starting indices
  integer n                     ! indices
  integer nlon_ice(plat)        ! number of lons per lat on bdy dataset
  integer i                     ! index into chunk
  integer j                     ! latitude index
  integer k
  integer ncol
  integer istat                 ! error return
  integer lchnk           ! chunk to process
  integer  :: yr, mon, day      ! components of a date
  integer  :: ncdate            ! current date in integer format [yyyymmdd]
  integer  :: ncsec             ! current time of day [seconds]
  real(r8) calday               ! calendar day (includes yr if no cycling)
  real(r8) caldayloc            ! calendar day (includes yr if no cycling)
  real(r8) xvar(plon,plat,2)    ! work space 
!-----------------------------------------------------------------------
!
! For aqua_planet there is no ice anywhere
!
     if(aqua_planet)then
        do lchnk=begchunk,endchunk
           aice(:pcols,lchnk) = 0.0
           previcefrac(:pcols,lchnk)=0.0
           call update_ocnice (lchnk)
        end do
        return
     end if
!
! Initialize time indices
!
  nm = 1
  np = 2
!
! Allocate space for data.
!
  allocate( ice(pcols,begchunk:endchunk), stat=istat )
  call alloc_err( istat, 'iceini', 'ice', &
       pcols*(endchunk-begchunk+1) )


  allocate( icebdy(pcols,begchunk:endchunk,2), stat=istat )
  call alloc_err( istat, 'iceini', 'icebdy', &
       pcols*(endchunk-begchunk+1)*2 )
!
! SPMD: Master does all the work.
!

  if (masterproc) then
!
! Use year information only if not cycling ice dataset
!
     if (is_first_step()) then
        dtime = get_step_size()
        dtime = -dtime
        calday = get_curr_calday(offset=dtime)
        call get_curr_date(yr, mon, day, ncsec,offset=dtime)
     else
        calday = get_curr_calday()
        call get_curr_date(yr, mon, day, ncsec)
     endif
     if (icecyc) then
        caldayloc = calday
     else
        caldayloc = calday + yr*daysperyear
     end if
     ncdate = yr*10000 + mon*100 + day
!
! Get and check dimension info
!
     call wrap_inq_dimid( ncid_sst, 'lon', londimid   )
     call wrap_inq_dimid( ncid_sst, 'time', timeid  )
     call wrap_inq_dimid( ncid_sst, 'lat', latdimid   )

     call wrap_inq_dimlen( ncid_sst, londimid, lonsiz   )
     if (lonsiz /= plon) then
        write(6,*)'ICEINI: lonsiz=',lonsiz,' must = plon=',plon
        call endrun
     end if
     call wrap_inq_dimlen( ncid_sst, latdimid, latsiz   )
     if (latsiz /= plat) then
        write(6,*)'ICEINI: latsiz=',latsiz,' must = plat=',plat
        call endrun
     end if
     call wrap_inq_dimlen( ncid_sst, timeid, timesiz   )
!
! Check to make sure space allocated for time variables is sufficient
!
     if (timesiz>toticesz) then
        write(6,*)'ICEINI:  Allocated space for ice data is insufficient.'
        write(6,*)'Please increase parameter toticesz to',timesiz,' and recompile.'
        call endrun
     end if
!
! Check to ensure reduced or not grid of dataset matches that of model
!
     if (fullgrid) then
        call wrap_inq_varid( ncid_sst, 'lon', lonid   )
     else
        call wrap_inq_varid (ncid_sst, 'nlon', nlonid)
        call wrap_get_var_int (ncid_sst, nlonid, nlon_ice)
        do j=1,plat
           if (nlon_ice(j) /= nlon(j)) then
              write(6,*)'ICEINI: model grid does not match dataset grid'
              call endrun
           end if
        end do
     end if

     call wrap_inq_varid( ncid_sst, 'date', dateid   )
     call wrap_inq_varid( ncid_sst, 'datesec', secid   )
     call wrap_inq_varid( ncid_sst, 'ice_cov', iceid   )
     call wrap_inq_varid( ncid_sst, 'lat', latid   )
!
! Retrieve entire date and sec variables.
!
     call wrap_get_var_int (ncid_sst,dateid,date_ice)
     call wrap_get_var_int (ncid_sst,secid,sec_ice)
     if (icecyc) then
        if (timesiz<12) then 
           write(6,*)'ICEINI: ERROR' 
           write(6,*)'When cycling ice, ice data set must have 12' 
           write(6,*)'consecutive months of data starting with Jan'
           write(6,*)'Current dataset has only ',timesiz,' months'
           call endrun
        end if
        do n = 1,12
           if (mod(date_ice(n),10000)/100/=n) then
              write(6,*)'ICEINI: ERROR' 
              write(6,*)'When cycling ice, ice data set must have 12' 
              write(6,*)'consecutive months of data starting with Jan'
              write(6,*)'Month ',n,' of ice data set is out of order'
              call endrun
           end if
        end do
     end if

     strt3(1) = 1
     strt3(2) = 1
     strt3(3) = 1
     cnt3(1)  = lonsiz
     cnt3(2)  = latsiz
     cnt3(3)  = 1
!
! Special code for interpolation between December and January
!
     if (icecyc) then
        n = 12
        np1 = 1
        call bnddyi(date_ice(n  ), sec_ice(n  ), cdayicem)
        call bnddyi(date_ice(np1), sec_ice(np1), cdayicep)
        if (caldayloc<=cdayicep .or. caldayloc>cdayicem) then
           strt3(3) = n
           call wrap_get_vara_realx (ncid_sst,iceid,strt3,cnt3,xvar(1,1,nm))
           strt3(3) = np1                                      
           call wrap_get_vara_realx (ncid_sst,iceid,strt3,cnt3,xvar(1,1,np))
           goto 10
        end if
     end if
!
! Normal interpolation between consecutive time slices.
!
     do n=1,timesiz-1
        np1 = n + 1
        call bnddyi(date_ice(n  ), sec_ice(n  ), cdayicem)
        call bnddyi(date_ice(np1), sec_ice(np1), cdayicep)
        if (.not.icecyc) then
           yr = date_ice(n)/10000
           cdayicem = cdayicem + yr*daysperyear
           yr = date_ice(np1)/10000
           cdayicep = cdayicep + yr*daysperyear
        end if
        if (caldayloc>cdayicem .and. caldayloc<=cdayicep) then
           strt3(3) = n
           call wrap_get_vara_realx (ncid_sst,iceid,strt3,cnt3,xvar(1,1,nm))
           strt3(3) = np1                                      
           call wrap_get_vara_realx (ncid_sst,iceid,strt3,cnt3,xvar(1,1,np))
           goto 10
        end if
     end do
     write(6,*)'ICEINI: Failed to find dates bracketing ncdate, ncsec=', ncdate, ncsec
     call endrun
10   continue
     write(6,*)'ICEINI: Read ice data for dates ',date_ice(n),sec_ice(n), &
          ' and ',date_ice(np1),sec_ice(np1)

#if (defined SPMD )
     call mpibcast( timesiz, 1, mpiint, 0, mpicom )
     call mpibcast( date_ice, toticesz, mpiint, 0, mpicom )
     call mpibcast( sec_ice, toticesz, mpiint, 0, mpicom )
     call mpibcast( cdayicem, 1, mpir8, 0, mpicom )
     call mpibcast( cdayicep, 1, mpir8, 0, mpicom )
     call mpibcast( np1, 1, mpiint, 0, mpicom )
  else
     call mpibcast( timesiz, 1, mpiint, 0, mpicom )
     call mpibcast( date_ice, toticesz, mpiint, 0, mpicom )
     call mpibcast( sec_ice, toticesz, mpiint, 0, mpicom )
     call mpibcast( cdayicem, 1, mpir8, 0, mpicom )
     call mpibcast( cdayicep, 1, mpir8, 0, mpicom )
     call mpibcast( np1, 1, mpiint, 0, mpicom )
#endif
  end if

  call scatter_field_to_chunk(1,1,2,plon,xvar,icebdy)
!
! Uncomment the following to get fractional model to give bit for bit 
! with non fractional
!
!!$  do lchnk=begchunk,endchunk
!!$     ncol = get_ncols_p(lchnk)
!!$     do i=1,ncol
!!$        do k=1,2
!!$           if (icebdy(i,lchnk,k)<=tsice) then
!!$              icebdy(i,lchnk,k) = 1.0
!!$           else
!!$              icebdy(i,lchnk,k) = 0.
!!$           end if
!!$        end do
!!$     end do
!!$  end do

  return
end subroutine iceini

!----------------------------------------------------------------------- 
! 
! BOP
!
! !IROUTINE: iceint
!
! !DESCRIPTION:
!
! if "aqua_planet", specify ICE's analytically (Jerry Olson).
! Otherwise, time interpolate ICE's to current time, reading in new monthly data if
! necessary.
!
! Method: 
! 
! Author: L.Bath
! 
!-----------------------------------------------------------------------
!
! !INTERFACE:
!
subroutine iceint(prev_timestep, startup)
!
! !USES:
!
  use rgrid, only: nlon
  use comsrf, only: snowhice, sicthk,surface_state2d,landfrac
  use time_manager, only: is_first_step, get_curr_date, get_curr_calday, get_step_size

#if ( defined SPMD )
  use mpishorthand, only: mpicom, mpir8
#endif
!
! !PARAMETERS:
!
  logical, intent(in), optional :: prev_timestep ! If using previous timestep, set to true
  logical, intent(in), optional :: startup       ! flag says to initialize previcefrac
!
! EOP
!
!---------------------------Common blocks-------------------------------
#include <comctl.h>
#include <comlun.h>
!---------------------------Local variables-----------------------------
  integer cnt3(3)        ! array of counts for each dimension
  integer strt3(3)       ! array of starting indices
  integer i,j,lchnk      ! indices
  integer ncol           ! number of columns in current chunk
  integer ntmp           ! temporary
  integer :: dtime       ! timestep size [seconds]
  real(r8) fact1, fact2  ! time interpolation factors
  integer :: yr, mon, day! components of a date
  integer :: ncdate      ! current date in integer format [yyyymmdd]
  integer :: ncsec       ! current time of day [seconds]
  real(r8) :: calday     ! current calendar day
  real(r8) caldayloc     ! calendar day (includes yr if no cycling)
  real(r8) deltat        ! time (days) between interpolating ice data
!
! Aqua planet variables
!
  real(r8) pi            ! 3.14159...
  real(r8) pio180        ! pi/180.
  real(r8) tmp           ! temporary
  real(r8) tmp1          ! temporary
  real(r8) t0_max        ! max reference temperature
  real(r8) t0_min        ! min reference temperature
  real(r8) t0_max6       ! max asymmetric reference temperature for option 6
  real(r8) t0_max7       ! max asymmetric reference temperature for option 7
  real(r8) maxlat        ! cutoff latitude poleward of which ICE = 0 deg C
  real(r8) shift         ! number of degrees peak ICE is shifted off equator
  real(r8) shift9        ! number of degrees peak ICE is shifted off equator for opt. 9
  real(r8) shift10       ! number of degrees peak ICE is shifted off equator for opt. 10
  real(r8) latcen        ! center of asymmetric ICE forcing
  real(r8) latrad6       ! radius of asymmetric ICE forcing for option 6
  real(r8) latrad8       ! radius of asymmetric ICE forcing for option 8
  real(r8) loncen        ! center of asymmetric ICE forcing
  real(r8) lonrad        ! radius of asymmetric ICE forcing
  real(r8) xvar(plon,plat,2)    ! work space 
  integer  ice_option    ! option of analytical ICE algorithm
  logical :: previous              ! If using previous timestep, set to true
  logical :: startup_lcl           ! startup
!
!-----------------------------------------------------------------------
     if(aqua_planet) return
!
! SPMD: Master does all the work.  Sends needed info to slaves
!
!
! Use year information only if a multiyear dataset
!
     if ( .not. present(prev_timestep) ) then
        previous = .false.
     else
        previous = prev_timestep
     end if

     if ( .not. present(startup) ) then
        startup_lcl = .false.
     else
        startup_lcl = startup
     end if
     
     if (previous .and. is_first_step()) then
        dtime = get_step_size()
        dtime = -dtime
        calday = get_curr_calday(offset=dtime)
        call get_curr_date(yr, mon, day, ncsec,offset=dtime)
     else
        calday = get_curr_calday()
        call get_curr_date(yr, mon, day, ncsec)
     endif
     if (icecyc) then
        caldayloc = calday
     else
        caldayloc = calday + yr*daysperyear
     end if

     if (masterproc) then

        strt3(1) = 1
        strt3(2) = 1
        strt3(3) = 1
        cnt3(1)  = lonsiz
        cnt3(2)  = latsiz
        cnt3(3)  = 1

     endif
!
! If model time is past current forward ice timeslice, read in the next
! timeslice for time interpolation.  Messy logic is for icecyc = .true. 
! interpolation between December and January (np1==1).  Note that 
! np1 is never 1 when icecyc is .false.
!
     if (caldayloc > cdayicep .and. .not. (np1==1 .and. caldayloc>cdayicem)) then
        if (icecyc) then
           np1 = mod(np1,12) + 1
        else
           np1 = np1 + 1
        end if
        if (np1>timesiz) then
           if (masterproc) then
              write(6,*)'ICEINT: Attempt to read past end of ICE dataset'
           endif
           call endrun
        end if
        cdayicem = cdayicep
        call bnddyi(date_ice(np1), sec_ice(np1), cdayicep)

        if (.not.icecyc) then
           yr = date_ice(np1)/10000
           cdayicep = cdayicep + yr*daysperyear
        end if

        if (np1==1 .or. caldayloc<=cdayicep) then
           ntmp = nm
           nm = np
           np = ntmp
           if (masterproc) then
              strt3(3) = np1
              call wrap_get_vara_realx (ncid_sst,iceid,strt3,cnt3,xvar(1,1,np))
              write(6,*)'ICEINT: Read ice for date (yyyymmdd) ',date_ice(np1), &
                   ' sec ',sec_ice(np1)
           endif
           call scatter_field_to_chunk(1,1,1,plon,xvar(1,1,np),icebdy(1,begchunk,np))
        else
           if (masterproc) then
              write(6,*)'ICEINT: Input ice for date',date_ice(np1), &
                   ' sec ',sec_ice(np1), 'does not exceed model date',ncdate,&
                   ' sec ',ncsec,' Stopping.'
           endif
           call endrun
        end if
     end if
!
! Time interpolation.  Account for December-January interpolation if
! cycling ice dataset.  Again note that np1 is never 1 when icecyc is false
!
     if (np1==1) then                    ! Dec-Jan interpolation
        deltat = cdayicep + daysperyear - cdayicem
        if (caldayloc>cdayicep) then      ! We're in December
           fact1 = (cdayicep + daysperyear - caldayloc)/deltat
           fact2 = (caldayloc - cdayicem)/deltat
        else                                ! We're in January
           fact1 = (cdayicep - caldayloc)/deltat
           fact2 = (caldayloc + daysperyear - cdayicem)/deltat
        end if
     else
        deltat = cdayicep - cdayicem
        fact1 = (cdayicep - caldayloc)/deltat
        fact2 = (caldayloc - cdayicem)/deltat
     end if
!
! Check sanity of time interpolation calculation to within 32-bit roundoff
!
     if (abs(fact1+fact2-1.)>1.e-6 .or. &
             fact1>1.000001 .or. fact1<-1.e-6 .or. &
             fact2>1.000001 .or. fact2<-1.e-6) then 
        if (masterproc) then
           write(6,*)'ICEINT: Bad fact1 and/or fact2=',fact1,fact2
        endif
        call endrun
     end if

     do lchnk=begchunk,endchunk
        ncol = get_ncols_p(lchnk)
        do i=1,ncol
              aice(i,lchnk) = icebdy(i,lchnk,nm)*fact1 + icebdy(i,lchnk,np)*fact2
              aice(i,lchnk)=min(aice(i,lchnk),1.)
              aice(i,lchnk)=max(aice(i,lchnk),0.)

              if (landfrac(i,lchnk) >= 1.) then
                aice(i,lchnk) = 0.
              end if
        end do
!
! get new ocean/ice percentages since ice coverage was just computed for this timestep
!
        call update_ocnice(lchnk)
!
! On startup, set previcefrac such that new ice will not be flagged in newiceproperties
!
        if (startup_lcl) then
           do i=1,ncol
              previcefrac(i,lchnk) = icefrac(i,lchnk)
           end do
        end if
     end do
!
! check for new ice and give it some default values for sicthk and snowhice
! and tssub
!
     call newiceproperties ()

     return
   end subroutine iceint

!======================================================================
! PRIVATE ROUTINES: Following routines are Private to this module
!======================================================================
!----------------------------------------------------------------------- 
! 
! BOP
!
! !IROUTINE: newiceproperties
!
! Purpose: 
!	Set the initial sea-ice surface variables associated with ICE.
!	This routine is private to this module.
!
!-----------------------------------------------------------------------
!
! !INTERFACE:
!
subroutine newiceproperties
!
! !USES:

  use comsrf, only: snowhice, sicthk, surface_state2d
  use sst_data, only: sst
  use time_manager, only: is_first_step
!
! EOP
!
!---------------------------Local variables-----------------------------
  integer lchnk, i,k  ! loop indices
  integer ncol        ! number of columns in current chunk
  real(r8) tsurf
!
! Set initial ice surface values
!
  do lchnk=begchunk,endchunk
     ncol = get_ncols_p(lchnk)
     do i=1,ncol
        if (icefrac(i,lchnk) > 0._r8) then
           snowhice(i,lchnk) = 0.005_r8
           sicthk(i,lchnk) = 2.0_r8
           if ( previcefrac(i,lchnk) == 0.) then !newice
              do k=1,plevmx
                 surface_state2d(lchnk)%tssub(i,k)=tsice+tmelt
              end do
           end if
        else !tssub for non ice boxes filled with open ocean temp
           snowhice(i,lchnk) = 0._r8
           sicthk(i,lchnk) = 0._r8
           do k=1,plevmx
              surface_state2d(lchnk)%tssub(i,k)=sst(i,lchnk)+tmelt
           end do
        end if
     end do
  end do
  return
end subroutine newiceproperties
end module ice_data
