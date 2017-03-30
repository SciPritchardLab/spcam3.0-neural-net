subroutine driver (ncprec)
!------------------------------------------------------------------
! Purpose:
!   Horizontally interpolate file containing aerosol masses to CAM grid.
!   Sum the values in the vertical such that output arrays contain the
!   column sum from the bottom of the atmosphere to that level.
!   Convert monthly averages to mid month values.
!   Input file assumed to be MATCH ncep runs, averaged by month.
!        and backsolved to provide Mid-month values)
!     
!  Method:
!    read data from file
!    interpolate data onto CAM horizontal grid
!
!------------------------------------------------------------------
   use prec
   use globals
   use preserve_mean, only: monthly_to_midmonth
   
   implicit none

   include 'netcdf.inc'
!
! Arguments
!
   integer, intent(in) :: ncprec     ! specify 32-bit or 64-bit precision for output variables
!
! Local workspace
!
   character(len=8), parameter :: aerosol_name(naer) =  &
        (/"MSUL    "&
         ,"MSSLT   "&
         ,"MDUST1  "&
         ,"MDUST2  "&
         ,"MDUST3  "&
         ,"MDUST4  "&
         ,"MOCPHO  "&
         ,"MBCPHO  "&
         ,"MOCPHI  "&
         ,"MBCPHI  "/)
   
   integer, parameter :: idxsslt = 2   ! index of SSLT
   
   integer :: lonidi = -1              ! longitude id input file
   integer :: latidi = -1              ! latitude id input file
   integer :: dateidi = -1             ! date id input file
   integer :: datesecidi = -1          ! datesec id input file
   integer :: mhybiid = -1             ! MATCH hybi id input file
   integer :: mpsid = -1               ! MATCH PS id input file
   integer :: species_id(naer)         ! aerosol ids input file

   integer :: lonido = -1              ! longitude id output file
   integer :: latido = -1              ! latitude id output file
   integer :: dateido = -1             ! date id output file
   integer :: datesecido = -1          ! datesec id output file
   integer :: hybiido = -1             ! MATCH hybi id output file
   integer :: psido = -1               ! MATCH PS id output file
   integer :: varido(naer)             ! aerosol ids output file
   integer :: timeido                  ! time id output filep

   integer :: nlonid                   ! nlon id output file (if present)
   integer :: rlonid                   ! rlon id output file (if present)

   real(r8) :: tempin(nxi,nz,nyi)      ! temp variable for bilinear interpolation
   real(r8) :: lonin(nxi)              ! longitude on input grid (rectangular)
   real(r8) :: latin(nyi)              ! latitude on input grid
   real(r8) :: m_hybi(nz+1)            ! MATCH hybi
   real(r8) :: m_ps(nxi,nyi,ntime)     ! surface pressure from MATCH

   real(r8) :: tempout(nxo,nz,nyo)     ! temp variable for bilinear interpolation
   real(r8) :: rlonout(nxo,nyo)        ! longitude on output grid (2-d for reduced grid)
   real(r8) :: lonout(nxo)             ! longitude on output grid (if rectangular)
   real(r8) :: latout(nyo)             ! latitude on input grid
   real(r8) :: M_ps_cam(nxo,nyo,ntime) ! surface pressure from MATCH on cam grid

   integer :: i, j, k                  ! x,y,z indices
   integer :: m                        ! constituent index
   integer :: mo                       ! month index
   integer :: dimids(nf_max_var_dims)  ! variable shape on output netcdf file
   integer :: date(ntime)              ! date (yyyymmdd)
   integer :: datesec(ntime)           ! seconds of date
   integer :: nlonout(nyo)             ! number of longitudes per latitude
   integer :: start(4)                 ! starting position in netcdf file
   integer :: kount(4)                 ! number of values to take from netcdf file
!
! a temporary place to store mmr's from files (generated from MATCH runs)
!
   real(r8) :: fspecies(nxi,nyi,nz)          ! aerosol mmr's from MATCH file
   real(r8), allocatable :: aerosol(:,:,:,:) ! aerosol mmr's from MATCH file on CAM grid
                                             ! allocate instead of dimension directly to prevent stack overflow
   allocate (aerosol(nxo,nyo,nz,ntime))
!
! Get required info from input file.
!
   call wrap_inq_varid (ncidi, 'lon', lonidi)
   call wrap_inq_varid (ncidi, 'lat', latidi)
   call wrap_inq_varid (ncidi, 'hybi', mhybiid)
   call wrap_inq_varid (ncidi, 'PS', mpsid)
   call wrap_inq_varid (ncidi, 'date', dateidi)
   call wrap_inq_varid (ncidi, 'datesec', datesecidi)

   call wrap_get_var_double (ncidi, lonidi, lonin)
   call wrap_get_var_double (ncidi, latidi, latin)
   call wrap_get_var_double (ncidi, mhybiid, m_hybi)
   call wrap_get_var_double (ncidi, mpsid, m_ps)
   call wrap_get_var_int (ncidi, dateidi, date)
   call wrap_get_var_int (ncidi, datesecidi, datesec)
!
! If time variable is not found on output file, create it.
! Then create date and datesec variables.
!
   if (nf_inq_varid (ncido, 'time', timeido) /= nf_noerr) then
      call wrap_def_var (ncido, 'time', nf_double, 1, timedimido, timeido)
   end if
   call wrap_def_var (ncido, 'date', nf_int, 1, timedimido, dateido)
   call wrap_def_var (ncido, 'datesec', nf_int, 1, timedimido, datesecido)
!
! Define hybi and PS on output grid
!
   call wrap_def_var (ncido, 'hybi', nf_double, 1, ilevdimido, hybiido)
   dimids(1:4) = (/londimido, latdimido, timedimido, -1/)
   call wrap_def_var (ncido, 'PS', ncprec, 3, dimids, psido)
!
! Read input variable names and define output names accordingly
! Append '_V' to indicate field has been vertically summed from sfc to each level
!
   dimids(1:5) = (/londimido, latdimido, levdimido, timedimido, -1/)
   do m = 1, naer
      call wrap_inq_varid (ncidi, trim (aerosol_name(m)), species_id(m))
      call wrap_def_var (ncido, trim (aerosol_name(m))//'_V', ncprec, 4, dimids, varido(m))
   end do
!
! Define global attribute "cam-ready" which will be checked for by CAM to prevent
! the use of datasets not run through the interpaerosols procedure.
!
   call wrap_put_att_text (ncido, nf_global, 'cam-ready', 'yes')
!
! End define mode on output file.  Copy required data from input file to output file
!
   call wrap_enddef (ncido)

   call wrap_put_var_double (ncido, hybiido, m_hybi)
!
! Retrieve output grid definition
!
   call wrap_inq_varid (ncido, 'lat', latido)
   call wrap_get_var_double (ncido, latido, latout)
!
! Reduced grid possibility: If coordinate variable 'lon' is on the dataset then
! the grid is rectangular.  If not, assume the grid is reduced and retrieve
! 2-d variable 'rlon' which is used to define the reduced grid.
!
   if (nf_inq_varid (ncido, 'lon', lonido) == nf_noerr) then
      call wrap_get_var_double (ncido, lonido, lonout)
      nlonout(1:nyo) = nxo
      do j=1,nyo
         rlonout(:,j) = lonout(:)
      end do
   else
      call wrap_inq_varid (ncido, 'nlon', nlonid)
      call wrap_inq_varid (ncido, 'rlon', rlonid)
      call wrap_get_var_int (ncido, nlonid, nlonout)
      call wrap_get_var_double (ncido, rlonid, rlonout)
   end if
!
! interpolate match's surface pressure and get mid-month values
!
   do mo=1,ntime
      call bilin (M_ps(1,1,mo), lonin, latin, nxi, nxi, &
                  1, 1, nyi, M_ps_cam(1,1,mo), rlonout, &
                  latout, nxo, nlonout, 1, nyo)
   end do
   call monthly_to_midmonth (M_ps_cam, 1)
!
! Retrieve Aerosol Masses (kg/m^2 in each layer)
!
   do m=1,naer
      do mo=1,ntime
         start(:) = (/1,1,1,mo/)
         kount(:) = (/nxi,nyi,nz,1/)
         call wrap_get_vara_double (ncidi, species_id(m), start, kount, fspecies)
!
! Accumulate mass below each interface level
!  note that lowest level (nz+1) is assumed to be zero 
!  but there isn't even storage for this value in the array
!
         do k=nz-1,1,-1
            fspecies(:,:,k) = fspecies(:,:,k) + fspecies(:,:,k+1)
         end do
!
! Transpose coords of retrieved aerosols to enable using CAM's bi-linear
! interpolation code.  Interpolate onto CAM horizontal grid
!
         do k=1,nz
            tempin(:,k,:) = fspecies(:,:,k)
         end do

         call bilin (tempin, lonin, latin, nxi, nxi, &
                     nz, nz, nyi, tempout, rlonout, &
                     latout, nxo, nlonout, nz, nyo)

         do k=1,nz
            aerosol(:,:,k,mo) = tempout(:,k,:)
         end do
!
! Sea Salt over land is minuscule.  After interpolation, differencing
! different levels of the cumulative mass can lead to underflow errors.
! To solve this problem, set sea salt to 0 for any column where
! the total mass is less than a threashold which looks like roundoff
! (and is unmeasurable?)
!
         if (m == idxSSLT) then
            do j=1,nyo
               do i=1,nxo
                  if (aerosol(i,j,1,mo) < 1.e-24) then
                     aerosol(i,j,:,mo) = 0.
                  end if
               end do
            end do
         end if
      end do ! mo
!
! convert from monthly average to mid-month values
!
      call monthly_to_midmonth (aerosol, nz)

      do mo=1,ntime
         do j=1,nyo
!
! make sure total column mass total is not negative
!
            do i=1,nxo
               aerosol(i,j,1,mo) = max (aerosol(i,j,1,mo), 0._r8)
            end do
!
! make function non-increasing and positive
!
            do k=2,nz
               do i=1,nxo
                  aerosol(i,j,k,mo) = min (aerosol(i,j,k,mo), aerosol(i,j,k-1,mo))
                  aerosol(i,j,k,mo) = max (aerosol(i,j,k,mo), 0._r8)
               end do
            end do
         end do
!
! Write interpolated data to output file
!
         start(:) = (/1,1,1,mo/)
         kount(:) = (/nxo,nyo,nz,1/)
         call wrap_put_vara_double (ncido, varido(m), start, kount, aerosol(1,1,1,mo))
      end do     ! loop over months (mo)
   end do        ! loop over constituents (m)

   do mo=1,ntime
      start(:) = (/1,1,mo,-1/)
      kount(:) = (/nxo,nyo,1,-1/)
      call wrap_put_vara_double (ncido, psido, start, kount, M_ps_cam(1,1,mo))

      start(:) = (/mo,-1,-1,-1/)
      kount(:) = (/1,-1,-1,-1/)
      call wrap_put_vara_int (ncido, dateido, start, kount, date(mo))
      call wrap_put_vara_int (ncido, datesecido, start, kount, datesec(mo))
   end do

   deallocate (aerosol)
   return
end subroutine driver
