#include <misc.h>
#include <preproc.h>

module histFileDGVMMod

#if (defined DGVM)

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: histFileDGVM
! 
! !DESCRIPTION: 
! Methods for DGVM history file handling
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varcon, only : spval  ! fill value
  use fileutils, only : get_filename
!
! !PUBLIC TYPES:
  implicit none
  save
  character(len=256), public :: dgvm_fn   !dgvm history filename
  integer , public :: ncid                !netcdf file id
  integer , public :: beg3d(3)            !netCDF 3-d start index 
  integer , public :: len3d(3)            !netCDF 3-d count index 
  integer , public :: beg4d(4)            !netCDF 4-d start index 
  integer , public :: len4d(4)            !netCDF 4-d count index 
!
! !PUBLIC MEMBER FUNCTIONS:
  public  :: histcrt_dgvm
  public  :: histwrt_dgvm
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: set_dgvm_filename            !sets dgvm history filename
!
! !REVISION HISTORY:
! Module created by Mariana Vertenstein
!
!EOP
!
! !PRIVATE TYPES:
  private
!
! Grid and time related variable id's
  integer :: mcdate_id        !id current date, yyyymmdd format (mcdate)
  integer :: mcsec_id         !id current seconds in day (mcsec)
  integer :: mdcur_id         !id current day (from base day) (mdcur)
  integer :: mscur_id         !id current seconds of current day (mdcur)
  integer :: nstep_id         !id current nstep 
  integer :: timcom_id        !id time comment (timcom)
  integer :: lonvar_id        !id full grid longitude coordinate variable
  integer :: latvar_id        !id full grid latitude  coordinate variable
  integer :: timvar_id        !id timecoordinate variable
  integer :: longxy_id        !id 2d longitudes (longxy)
  integer :: latixy_id        !id 2d latitudes (latixy)
  integer :: numlon_id        !id number of longitudes at each latitude
  integer :: landmask_id      !id 2d land/ocean mask (landmask)
#if (defined OFFLINE)
  integer :: edgen_id         !id northern edge of grid (lsmedge(1))
  integer :: edgee_id         !id eastern  edge of grid (lsmedge(2))
  integer :: edges_id         !id southern edge of grid (lsmedge(3))
  integer :: edgew_id         !id western  edge of grid (lsmedge(4))
#endif
  integer :: ncprec           !output precision  
!
! Output variable id's
  integer, public :: afirefrac_id         !output variable id 
  integer, public :: acfluxfire_id        !output variable id  
  integer, public :: bmfm_id              !output variable id  
  integer, public :: afmicr_id            !output variable id 
  integer, public :: fpcgrid_id           !output variable id 
  integer, public :: itypveg_id           !output variable id  
  integer, public :: lmind_id             !output variable id 
  integer, public :: rmind_id             !output variable id 
  integer, public :: smind_id             !output variable id 
  integer, public :: hmind_id             !output variable id 
  integer, public :: nind_id              !output variable id 
  integer, public :: begwater_id          !output variable id  
  integer, public :: endwater_id          !output variable id 
  integer, public :: begenergy_id         !output variable id    
  integer, public :: endenergy_id         !output variable id  
!----------------------------------------------------------------------- 

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: histcrt_dgvm
!
! !INTERFACE:
  subroutine histcrt_dgvm ()
!
! !DESCRIPTION: 
! create netCDF history file
! This subroutine opens a new netCDF dgvm history file. 
! Global attributes and variables are defined in define mode. 
! Upon exiting this routine, define mode is exited and the file 
! is ready to write.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clm_varpar, only : lsmlon, lsmlat, maxpatch_pft
    use clm_varctl, only : ctitle, finidat, fsurdat, fpftcon, &
         frivinp_rtm, caseid
    use clm_varsur, only : fullgrid, offline_rdgrid
    use time_manager, only : get_ref_date	
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
!
! !CALLED FROM:
! subroutine lpj in module DGVMMod
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer i                   !field do loop index
    integer status              !netCDF error status
    integer lon_id              !netCDF id for longitude dimension
    integer lat_id              !netCDF id for latitude dimension
    integer pft_id              !netCDF id for patch dimension
    integer tim_id              !netCDF id for time dimension
    integer strlen_id           !netCDF id for character string variables
    integer dim1_id(1)          !netCDF dimension id for 1-d variables
    integer dim2_id(2)          !netCDF dimension id for 2-d variables
    integer dim3_id(3)          !netCDF dimension id for 3-d variables
    integer dim4_id(4)          !netCDF dimension id for 4-d variables
    integer omode               !netCDF dummy variable
    integer ret                 !netCDF return status 
    character(len=256) name     !name of attribute
    character(len=256) str      !temporary string
    character(len=  8) curdate  !current date
    character(len=  8) curtime  !current time 
    character(len= 10) basedate !base date (yyyymmdd)
    character(len=  8) basesec  !base seconds
    integer yr,mon,day,nbsec    !year,month,day,seconds components of a date
    integer hours,minutes,secs  !hours,minutes,seconds of hh:mm:ss
!-----------------------------------------------------------------------

    ! Set output precision

    ncprec = nf_double

    ! Create new netCDF file. File will be in define mode

    dgvm_fn = set_dgvm_filename()
    call wrap_create (trim(dgvm_fn), nf_clobber, ncid)

    ! Set fill mode to "no fill" to optimize performance

    status = nf_set_fill (ncid, nf_nofill, omode)
    if (status /= nf_noerr) then
       write(6,*) ' netCDF error = ',nf_strerror(status)
       call endrun
    end if

    ! Create global attributes. 

    str = 'CF1.0'
    call wrap_put_att_text (ncid, NF_GLOBAL, 'conventions', trim(str))

    call getdatetime(curdate, curtime)
    str = 'created on ' // curdate // ' ' // curtime
    call wrap_put_att_text(ncid, NF_GLOBAL,'history', trim(str))

    call getenv ('LOGNAME', str)
    call wrap_put_att_text (ncid, NF_GLOBAL, 'logname',trim(str))

    call getenv ('HOST', str)
    call wrap_put_att_text (ncid, NF_GLOBAL, 'host', trim(str))

    str = 'Community Land Model: CLM2'
    call wrap_put_att_text (ncid, NF_GLOBAL, 'source', trim(str))

    str = '$Name:  $'
    call wrap_put_att_text (ncid, NF_GLOBAL, 'version', trim(str))

    str = '$Id: histFileDGVMMod.F90,v 1.1.4.2 2003/06/13 16:13:44 hender Exp $'
    call wrap_put_att_text (ncid, NF_GLOBAL, 'revision_id', trim(str))

    str = ctitle 
    call wrap_put_att_text (ncid, NF_GLOBAL, 'case_title', trim(str))

    str = caseid
    call wrap_put_att_text (ncid, NF_GLOBAL, 'case_id', trim(str))

    str = 'created at run time'
    if (fsurdat /= ' ') str = get_filename(fsurdat)
    call wrap_put_att_text(ncid, NF_GLOBAL, 'Surface_dataset', trim(str))

    str = 'arbitrary initialization'
    if (finidat /= ' ') str = get_filename(finidat)
    call wrap_put_att_text(ncid, NF_GLOBAL, 'Initial_conditions_dataset', trim(str))

    str = get_filename(fpftcon)
    call wrap_put_att_text(ncid, NF_GLOBAL, 'PFT_physiological_constants_dataset', trim(str))

    if (frivinp_rtm /= ' ') then
       str = get_filename(frivinp_rtm)
       call wrap_put_att_text(ncid, NF_GLOBAL, 'RTM_input_datset', trim(str))
    endif

    ! Define dimensions. 

    call wrap_def_dim (ncid, 'lon', lsmlon, lon_id)
    call wrap_def_dim (ncid, 'lat', lsmlat, lat_id)
    call wrap_def_dim (ncid, 'pft', maxpatch_pft, pft_id)
    call wrap_def_dim (ncid, 'time', nf_unlimited, tim_id)
    call wrap_def_dim (ncid, 'string_length', 80, strlen_id)

    ! Define coordinate variables (including time)

    if (fullgrid) then
       dim1_id(1) = lon_id
       call wrap_def_var (ncid, 'lon' , nf_double, 1, dim1_id, lonvar_id)
       call wrap_put_att_text (ncid, lonvar_id, 'long_name','coordinate longitude')
       call wrap_put_att_text (ncid, lonvar_id, 'units'    ,'degrees_east')

       dim1_id(1) = lat_id
       call wrap_def_var (ncid, 'lat' , nf_double, 1, dim1_id, latvar_id)
       call wrap_put_att_text (ncid, latvar_id, 'long_name','coordinate latitude')
       call wrap_put_att_text (ncid, latvar_id, 'units'    ,'degrees_north')
    endif

    call get_ref_date(yr, mon, day, nbsec)
    hours   = nbsec / 3600
    minutes = (nbsec - hours*3600) / 60
    secs    = (nbsec - hours*3600 - minutes*60)
    write(basedate,80) yr,mon,day
80  format(i4.4,'-',i2.2,'-',i2.2)
    write(basesec ,90) hours, minutes, secs
90  format(i2.2,':',i2.2,':',i2.2)
    dim1_id(1) = tim_id
    call wrap_def_var (ncid, 'time', ncprec, 1, dim1_id, timvar_id)
    call wrap_put_att_text (ncid, timvar_id, 'long_name','time')
    str = 'days since ' // basedate // " " // basesec
    call wrap_put_att_text (ncid, timvar_id, 'units'    ,str)
    call wrap_put_att_text (ncid, timvar_id, 'calendar' ,'noleap')

#if (defined OFFLINE)

    ! Define surface grid edge variables

    if (.not. offline_rdgrid) then
       call wrap_def_var (ncid , 'edgen', ncprec, 0, 0, edgen_id)
       call wrap_put_att_text (ncid, edgen_id, 'long_name',&
            'northern edge of surface grid')
       call wrap_put_att_text (ncid, edgen_id, 'units',&
            'degrees_north')

       call wrap_def_var (ncid, 'edgee', ncprec,0, 0, edgee_id)
       call wrap_put_att_text (ncid, edgee_id, 'long_name',&
            'eastern edge of surface grid')
       call wrap_put_att_text (ncid, edgee_id, 'units',&
            'degrees_east')

       call wrap_def_var (ncid, 'edges', ncprec, 0, 0, edges_id)
       call wrap_put_att_text (ncid, edges_id, 'long_name',&
            'southern edge of surface grid')
       call wrap_put_att_text (ncid, edges_id, 'units',&
            'degrees_north')

       call wrap_def_var (ncid, 'edgew', ncprec, 0, 0, edgew_id)
       call wrap_put_att_text (ncid, edgew_id , 'long_name',&
            'western edge of surface grid')
       call wrap_put_att_text (ncid, edgew_id , 'units',&
            'degrees_east')
    endif

#endif

    ! Define longitude, latitude, surface type: real (lsmlon x lsmlat) variables

    dim2_id(1) = lon_id
    dim2_id(2) = lat_id

    if (fullgrid) then
       call wrap_def_var (ncid, 'longxy' , ncprec, 2, dim2_id, longxy_id)
       call wrap_put_att_text (ncid, longxy_id, 'long_name','longitude')
       call wrap_put_att_text (ncid, longxy_id, 'units'    ,'degrees_east')
    else
       call wrap_def_var (ncid, 'rlongxy', ncprec, 2, dim2_id, longxy_id)
       call wrap_put_att_text (ncid, longxy_id, 'long_name','rlongitude')
       call wrap_put_att_text (ncid, longxy_id, 'units'    ,'degrees_east')
    endif

    call wrap_def_var (ncid, 'latixy', ncprec, 2, dim2_id, latixy_id)
    call wrap_put_att_text (ncid, latixy_id, 'long_name','latitude')
    call wrap_put_att_text (ncid, latixy_id, 'units'    ,'degrees_north')

    ! Define number of longitudes per latitude (reduced grid only) variable

    dim1_id(1) = lat_id

    call wrap_def_var (ncid, 'numlon', nf_int, 1, dim1_id, numlon_id)
    call wrap_put_att_text (ncid, numlon_id, 'long_name', &
         'number of longitudes at each latitude')

    ! Number of surface type variable

    dim2_id(1) = lon_id
    dim2_id(2) = lat_id

    call wrap_def_var (ncid, 'landmask', nf_int, 2, dim2_id, landmask_id)
    call wrap_put_att_text (ncid,landmask_id,'long_name',&
         'land/ocean mask (0=ocean, 1=land)')

    ! Define time information variables

    dim1_id(1) = tim_id

    call wrap_def_var (ncid , 'mcdate', nf_int, 1, dim1_id  , mcdate_id)
    call wrap_put_att_text (ncid, mcdate_id, 'long_name',&
         'current date (YYYYMMDD)')

    call wrap_def_var (ncid , 'mcsec' , nf_int, 1, dim1_id , mcsec_id)
    call wrap_put_att_text (ncid, mcsec_id, 'long_name',&
         'current seconds of current date')
    call wrap_put_att_text (ncid, mcsec_id, 'units','s')

    call wrap_def_var (ncid , 'mdcur' , nf_int, 1, dim1_id , mdcur_id)
    call wrap_put_att_text (ncid, mdcur_id, 'long_name',&
         'current day (from base day)')

    call wrap_def_var (ncid , 'mscur' , nf_int, 1, dim1_id , mscur_id)
    call wrap_put_att_text (ncid, mscur_id, 'long_name',&
         'current seconds of current day')
    call wrap_put_att_text (ncid, mcsec_id, 'units','s')

    call wrap_def_var (ncid , 'nstep' , nf_int, 1, dim1_id , nstep_id)
    call wrap_put_att_text (ncid, nstep_id, 'long_name',&
         'time step')

    ! Define DGVM output variables

    dim3_id(1) = lon_id
    dim3_id(2) = lat_id
    dim3_id(3) = tim_id

    dim4_id(1) = lon_id
    dim4_id(2) = lat_id
    dim4_id(3) = pft_id
    dim4_id(4) = tim_id 

    !written from lpj  
    call wrap_def_var (ncid, 'BURN', ncprec, 3, dim3_id, afirefrac_id)
    call wrap_put_att_text (ncid, afirefrac_id, 'long_name' , &
         'fraction of vegetated area burned') 
    call wrap_put_att_text (ncid, afirefrac_id, 'units',&
         'fraction of vegetated area')
    call wrap_put_att_realx(ncid, afirefrac_id, '_FillValue', ncprec,1 ,spval)

    call wrap_def_var (ncid, 'CFLUXFIRE', ncprec, 3, dim3_id, acfluxfire_id)
    call wrap_put_att_text (ncid, acfluxfire_id, 'long_name' ,&
         'carbon flux to the atmosphere due to fire')
    call wrap_put_att_text (ncid, acfluxfire_id, 'units'     ,&
         'grams C/m2 of vegetated area')
    call wrap_put_att_realx(ncid, acfluxfire_id, '_FillValue', ncprec,1 ,spval)

    call wrap_def_var (ncid, 'NPP', ncprec, 4, dim4_id, bmfm_id)
    call wrap_put_att_text (ncid, bmfm_id, 'long_name' ,&
         'annual net primary production' )
    call wrap_put_att_text (ncid, bmfm_id, 'units'     ,&
         'grams C/m2 of patch')
    call wrap_put_att_realx(ncid, bmfm_id, '_FillValue', ncprec,1 ,spval)

    call wrap_def_var (ncid, 'Rh', ncprec, 4, dim4_id, afmicr_id)
    call wrap_put_att_text (ncid, afmicr_id, 'long_name' ,&
         'annual heterotrophic respiration')  
    call wrap_put_att_text (ncid, afmicr_id, 'units', &
         'grams C/m2 of vegetated area')
    call wrap_put_att_realx(ncid, afmicr_id, '_FillValue', ncprec,1 ,spval)

    !written from lpjreset2
    call wrap_def_var (ncid, 'PFT', nf_int, 4, dim4_id, itypveg_id)  
    call wrap_put_att_text (ncid, itypveg_id, 'long_name' ,&
         'plant functional type') 
    ret = nf_put_att_int (ncid, itypveg_id, '_FillValue', nf_int, 1, 9999)
    if (ret/=NF_NOERR) call handle_error (ret)

    call wrap_def_var (ncid, 'FPCGRID', ncprec, 4, dim4_id, fpcgrid_id)
    call wrap_put_att_text (ncid, fpcgrid_id, 'long_name' ,&
         'plant functional type cover')
    call wrap_put_att_text (ncid, fpcgrid_id, 'units',&
         'fraction of vegetated area')
    call wrap_put_att_realx(ncid, fpcgrid_id, '_FillValue', ncprec,1 ,spval)

    call wrap_def_var (ncid, 'LCIND', ncprec, 4, dim4_id, lmind_id)
    call wrap_put_att_text (ncid, lmind_id, 'long_name' ,&
         'leaf carbon per individual') 
    call wrap_put_att_text (ncid, lmind_id, 'units',&
         'grams C')
    call wrap_put_att_realx(ncid, lmind_id, '_FillValue', ncprec,1 ,spval)

    call wrap_def_var (ncid, 'RCIND', ncprec, 4, dim4_id, rmind_id)
    call wrap_put_att_text (ncid, rmind_id, 'long_name' ,&
         'root carbon per individual')
    call wrap_put_att_text (ncid, lmind_id, 'units',&
         'grams C')
    call wrap_put_att_realx(ncid, rmind_id, '_FillValue', ncprec,1 ,spval)

    call wrap_def_var (ncid, 'SCIND', ncprec, 4, dim4_id, smind_id)
    call wrap_put_att_text (ncid, smind_id, 'long_name',&
         'stem carbon per individual')
    call wrap_put_att_text (ncid, smind_id, 'units', &
         'grams C')
    call wrap_put_att_realx(ncid, smind_id, '_FillValue', ncprec,1 ,spval)

    call wrap_def_var (ncid, 'HCIND', ncprec, 4, dim4_id, hmind_id)
    call wrap_put_att_text (ncid, hmind_id, 'long_name' ,& 
         'heartwood carbon per individual') 
    call wrap_put_att_text (ncid, hmind_id, 'units', &
         'grams C')
    call wrap_put_att_realx(ncid, hmind_id, '_FillValue', ncprec,1 ,spval)

    call wrap_def_var (ncid, 'NIND', ncprec, 4, dim4_id, nind_id)
    call wrap_put_att_text (ncid, nind_id, 'long_name' ,&
         'number of individuals') 
    call wrap_put_att_text (ncid, nind_id, 'units',&
         'individuals/m2 vegetated land')
    call wrap_put_att_realx(ncid, nind_id, '_FillValue', ncprec,1 ,spval)

    call wrap_def_var (ncid, 'BEGWATER', ncprec, 3, dim3_id, begwater_id)
    call wrap_put_att_text (ncid, begwater_id, 'long_name' ,&
         'beginning water')  
    call wrap_put_att_text (ncid, begwater_id, 'units',&
         'kg/m2')
    call wrap_put_att_realx(ncid, begwater_id, '_FillValue', ncprec,1 ,spval)

    call wrap_def_var (ncid, 'ENDWATER', ncprec, 3, dim3_id, endwater_id)
    call wrap_put_att_text (ncid, endwater_id, 'long_name', &
         'ending water')  
    call wrap_put_att_text (ncid, endwater_id, 'units',&
         'kg/m2')
    call wrap_put_att_realx(ncid, endwater_id, '_FillValue', ncprec,1 ,spval)

    call wrap_def_var (ncid, 'BEGENERGY', ncprec, 3, dim3_id, begenergy_id)
    call wrap_put_att_text (ncid, begenergy_id, 'long_name',& 
         'beginning energy')  
    call wrap_put_att_text (ncid, begenergy_id, 'units',&
         'J/m2') 
    call wrap_put_att_realx(ncid, begenergy_id, '_FillValue', ncprec,1 ,spval)

    call wrap_def_var (ncid, 'ENDENERGY', ncprec, 3, dim3_id, endenergy_id)
    call wrap_put_att_text (ncid, endenergy_id, 'long_name',&
         'ending energy')  
    call wrap_put_att_text (ncid, endenergy_id, 'units',&
         'J/m2')
    call wrap_put_att_realx(ncid, endenergy_id, '_FillValue', ncprec,1 ,spval)

    ! Setup beginning and ending write indices

    beg3d(1) = 1  ; len3d(1) = lsmlon
    beg3d(2) = 1  ; len3d(2) = lsmlat
    beg3d(3) = 1  ; len3d(3) = 1

    beg4d(1) = 1  ; len4d(1) = lsmlon
    beg4d(2) = 1  ; len4d(2) = lsmlat
    beg4d(3) = 1  ; len4d(3) = maxpatch_pft
    beg4d(4) = 1  ; len4d(4) = 1

    ! Finish creating netCDF file (end define mode)

    status = nf_enddef(ncid)

  end subroutine histcrt_dgvm

  !=======================================================================

  subroutine histwrt_dgvm()

    use shr_kind_mod, only: r8 => shr_kind_r8	
    use clm_varpar, only : lsmlon, lsmlat	
    use clm_varsur, only : longxy, latixy, landmask, fullgrid, numlon, &
         lsmedge, offline_rdgrid
    use time_manager, only : get_nstep, get_curr_date, get_curr_time
    use spmdMod, only : masterproc	
    use shr_const_mod, only : SHR_CONST_CDAY
    implicit none
    include 'netcdf.inc'

    ! ------------------------ local variables ---------------------------
    integer  :: i,j,k,l,m,n                    !do loop indices
    integer  :: beg1d(1),len1d(1)              !netCDF 1-d start,count indices
    integer  :: beg2d(2),len2d(2)              !netCDF 2-d start,count indices
    real(r8) :: lonvar(lsmlon)                 !only used for full grid 
    real(r8) :: latvar(lsmlat)                 !only used for full grid
    real(r8) :: time                           !current time
    integer  :: mdcur, mscur, mcdate           !outputs from get_curr_time
    integer  :: yr,mon,day,mcsec               !outputs from get_curr_date
    integer  :: nstep                          !time step
    ! --------------------------------------------------------------------

    if (masterproc) then

       ! Write out time-invariant variables 

#if (defined OFFLINE)
       if (.not. offline_rdgrid) then
          call wrap_put_var_realx (ncid, edgen_id, lsmedge(1))
          call wrap_put_var_realx (ncid, edgee_id, lsmedge(2))
          call wrap_put_var_realx (ncid, edges_id, lsmedge(3))
          call wrap_put_var_realx (ncid, edgew_id, lsmedge(4))
       endif
#endif

       ! Surface grid (coordinate variables, latitude, longitude, surface type). 

       if (fullgrid) then
          lonvar(1:lsmlon) = longxy(1:lsmlon,1)
          call wrap_put_var_realx (ncid, lonvar_id, lonvar)
          latvar(1:lsmlat) = latixy(1,1:lsmlat)
          call wrap_put_var_realx (ncid, latvar_id, latvar)
       endif
       call wrap_put_var_realx (ncid, longxy_id  , longxy)
       call wrap_put_var_realx (ncid, latixy_id  , latixy)
       call wrap_put_var_int   (ncid, landmask_id, landmask)
       call wrap_put_var_int   (ncid, numlon_id  , numlon)

       ! Write current date, current seconds, current day, current nstep

       beg1d(1) = 1 ; len1d(1) = 1
       call get_curr_date(yr, mon, day, mcsec)
       mcdate = yr*10000 + mon*100 + day
       call get_curr_time(mdcur,mscur)  
       time = mdcur + mscur/SHR_CONST_CDAY
       nstep = get_nstep()

       call wrap_put_vara_int (ncid, mcdate_id, beg1d, len1d, mcdate)
       call wrap_put_vara_int (ncid, mcsec_id , beg1d, len1d, mcsec)
       call wrap_put_vara_int (ncid, mdcur_id , beg1d, len1d, mdcur)
       call wrap_put_vara_int (ncid, mscur_id , beg1d, len1d, mscur)
       call wrap_put_vara_int (ncid, nstep_id , beg1d, len1d, nstep)
       call wrap_put_vara_realx (ncid, timvar_id, beg1d, len1d, time)

    endif ! end of if-masterproc block 

  end subroutine histwrt_dgvm

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_dgvm_filename
!
! !INTERFACE:
  character(len=256) function set_dgvm_filename ()
!
! !DESCRIPTION: 
! Determine initial dataset filenames
!
! !USES:
    use clm_varctl  , only : caseid
    use time_manager, only : get_curr_date
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    character(len=256) :: cdate       !date char string
    integer :: day                    !day (1 -> 31)
    integer :: mon                    !month (1 -> 12)
    integer :: yr                     !year (0 -> ...)
    integer :: sec                    !seconds into current day
!-----------------------------------------------------------------------

    call get_curr_date (yr, mon, day, sec) 
    write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr,mon,day,sec
    set_dgvm_filename = "./"//trim(caseid)//".clm2.hv."//trim(cdate)//".nc"

  end function set_dgvm_filename

#endif

end module histFileDGVMMod




