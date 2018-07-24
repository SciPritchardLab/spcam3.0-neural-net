	subroutine crmhistory(lchnk, icol, time_crm, zs, prechist, qrs_crm, qrl_crm,
     &	                fsds, fsns, fsnt, fsut,
     &	                flds, flns, flut, fsntc, flutc, fsdsc, flnsc, ps)
        use shr_kind_mod, only: r8 => shr_kind_r8
        use phys_grid, only: get_rlon_p, get_rlat_p, get_lon_p, 
     &	                     get_lat_p
        use history,   only: get_mtapes,tape,ngroup, hfilename_spec,
     &                             sec2hms, date2yyyymmdd, ctitle
        use filenames, only: interpret_filename_spec
        use time_manager, only: get_step_size, get_ref_date,
     &                          get_curr_date,get_curr_time,get_prev_date
        use time_manager, only:	is_last_step
        use pmgrid, only: plon,plond, plat, plev,plevp, iam
        use commap, only: latdeg,londeg, clat,clon

	implicit none
        include 'netcdf.inc'
	include 'common.inc'

#include <comhyb.h>

	integer lchnk,icol
	 real, intent(in) :: time_crm ! seconds of CRM simulation this CAM step
         real(r8), intent(in) :: qrs_crm(crm_nx,crm_ny,crm_nz)
         real(r8), intent(in) :: qrl_crm(crm_nx,crm_ny,crm_nz)
         real(r8), intent(in) :: prechist(crm_nx,crm_ny)! Precipitation rate at Surface
         real(r8), intent(in) :: fsds(crm_nx,crm_ny)   ! Shortwave Downwelling Flux at Surface
         real(r8), intent(in) :: fsns(crm_nx,crm_ny)   ! Surface solar absorbed flux
         real(r8), intent(in) :: fsnt(crm_nx,crm_ny)   ! Net solar flux at top of atmosphere
         real(r8), intent(in) :: fsut(crm_nx,crm_ny)   ! Shortwave Upwelling Flux at TOA
         real(r8), intent(in) :: flds(crm_nx,crm_ny)   ! Surface longwave downward flux
         real(r8), intent(in) :: flns(crm_nx,crm_ny)   ! Net upward longwave flux at surface
         real(r8), intent(in) :: flut(crm_nx,crm_ny)   ! Outgoing lw flux at model top
         real(r8), intent(in) :: fsntc(crm_nx,crm_ny)  ! Clearsky net solar flux at TOA
         real(r8), intent(in) :: fsdsc(crm_nx,crm_ny)  ! Clearsky downward solar flux at surface
         real(r8), intent(in) :: flnsc(crm_nx,crm_ny)  ! Clearsky net upward longwave flux at srf
         real(r8), intent(in) :: flutc(crm_nx,crm_ny)  ! Clearsky outgoing lw flux at model top
         real(r8), intent(in) :: ps ! surface pressure (Pa)
         real(r8), intent(in) :: zs ! surface elevation (m)


        integer, parameter :: max_chars = 128 ! max chars for char variables
	logical group
	integer grouplon1,grouplon2,grouplat1,grouplat2
	integer latcol,loncol ! lat, lon index for grid cell
        integer, parameter :: ptapes = 6             ! max number of tapes
	integer mtapes
	integer it,g,i,k,j
	integer lon,lat
	integer group_tape,group_number,group_lon,group_lat
	integer maxgroups ! number of groups of columns
	integer dimen5g(5),dimen6g(6),dimen1(1),dimen2t(2)
	integer start(6),count(6)
        integer dimenchar(2)       ! character dimension ids
        integer chardim            ! character dimension id
        integer tbnddim            ! time_bnds dimension id
	parameter (maxgroups=20)
	integer nfid(maxgroups,ptapes)
	integer nfils(maxgroups,ptapes) ! time sample count
	save nfid,nfils
      integer grouplondim(maxgroups) ! longitude dimension id
      integer grouplatdim(maxgroups) ! latitude dimension id
      integer :: ncreal          ! netCDF real data type
      integer :: mdtid            ! var id for timestep
      integer :: ndbaseid(maxgroups,ptapes)         ! var id for base day
      integer :: nsbaseid(maxgroups,ptapes)         ! var id for base seconds of base day
      integer :: nbdateid(maxgroups,ptapes)         ! var id for base date
      integer :: nbsecid(maxgroups,ptapes)          ! var id for base seconds of base date
      integer :: ndcurid(maxgroups,ptapes)          ! var id for current day
      integer :: nscurid(maxgroups,ptapes)          ! var id for current seconds of current day
      integer :: dateid(maxgroups,ptapes)           ! var id for current date
      integer :: datesecid(maxgroups,ptapes)        ! var id for curent seconds of current date
      integer :: timeid(maxgroups,ptapes)           ! var id for time
      integer :: tbndid(maxgroups,ptapes)           ! var id for time bounds          !
      integer :: date_writtenid(maxgroups,ptapes)   ! var id for date time sample written
      integer :: time_writtenid(maxgroups,ptapes)   ! var id for time time sample written
      integer :: nscurf(maxgroups,ptapes)           ! First "current" second of day for each h-file
      integer :: ncsecf(maxgroups,ptapes)           ! First "current" second of date for each h-file
      save ndbaseid,nsbaseid,nbdateid,nbsecid,ndcurid,nscurid
      save dateid,datesecid,timeid,date_writtenid,time_writtenid,tbndid
      save nscurf,ncsecf
      character(len=8) :: cdate  ! system date
      character(len=8) :: ctime  ! system time
      integer :: ndbase = 0      ! days component of base time
      integer :: nsbase = 0      ! seconds component of base time
      integer :: nbdate          ! base date in yyyymmdd format
      integer :: nbsec           ! time of day component of base date [seconds]
      integer :: iyr, imon, iday ! year, month, day components of a date
      integer ncsec ! seconds of current day
! varid for fields
      integer varid_u(maxgroups,ptapes)
      integer varid_v(maxgroups,ptapes)
      integer varid_w(maxgroups,ptapes)
      integer varid_t(maxgroups,ptapes)
      integer varid_qv(maxgroups,ptapes)
      integer varid_qc(maxgroups,ptapes)
      integer varid_qi(maxgroups,ptapes)
      integer varid_qr(maxgroups,ptapes)
      integer varid_qs(maxgroups,ptapes)
      integer varid_qrs(maxgroups,ptapes)
      integer varid_qrl(maxgroups,ptapes)
      integer varid_prec(maxgroups,ptapes)
      integer varid_fsds(maxgroups,ptapes)
      integer varid_fsns(maxgroups,ptapes)
      integer varid_fsnt(maxgroups,ptapes)
      integer varid_fsut(maxgroups,ptapes)
      integer varid_fsdsc(maxgroups,ptapes)
      integer varid_fsntc(maxgroups,ptapes)
      integer varid_flds(maxgroups,ptapes)
      integer varid_flns(maxgroups,ptapes)
      integer varid_flut(maxgroups,ptapes)
      integer varid_flnsc(maxgroups,ptapes)
      integer varid_flutc(maxgroups,ptapes)
      integer ps_var(maxgroups,ptapes)
      save varid_u,varid_v,varid_w,varid_t,varid_qv,varid_qc,varid_qi
      save varid_qr,varid_qs,varid_qrs,varid_qrl,varid_prec
      save varid_fsds,varid_fsns,varid_fsnt,varid_fsut
      save varid_flds,varid_flns,varid_flut
      save varid_fsdsc,varid_fsntc,varid_flnsc,varid_flutc,ps_var
      real(r8) :: tdata(2)         ! time interval boundaries
      real(r8) :: beg_time(maxgroups,ptapes) ! time at beginning of an averaging interval
      save beg_time

      integer, parameter :: nlen = 256    ! Length of strings
      character(len=nlen) :: nhfil(maxgroups,ptapes)    ! Array of current file names
      character(len=nlen) :: string    ! Temporary character string
      integer lenstring
      save nhfil
      integer ntapes
      character(len=max_chars) str ! character temporary
      character(len=8)  :: logname             ! user name
      character(len=16) :: host                ! host name
      integer crm_x_var
      integer crm_y_var
      integer crm_z_var
      integer pres_var, ps0_var
      integer z_var, zs_var
      integer crm_x_dim
      integer crm_y_dim
      integer crm_z_dim
      real(r8) crm_x(crm_nx)
      real(r8) crm_y(crm_ny)
      real(r8) presref(crm_nz) ! reference pressure assuming 1000 mb surface pressure
      real*4 work3d(crm_nx,crm_ny,crm_nz)
      real*4 work2d(crm_nx,crm_ny)
      integer :: ncdate            ! current date in integer format [yyyymmdd]
      integer :: ndcur             ! day component of current time
      integer :: nscur             ! seconds component of current time
      real*8 :: htime              ! current time
!
! netcdf variables
!
      integer ret                ! function return value
      integer timdim             ! unlimited dimension id
      integer latvar             ! latitude variable id
      integer lonvar             ! longitude variable id
      integer glatvar(maxgroups) ! group latitude variable id
      integer glonvar(maxgroups) ! group longitude variable id

      latcol = get_lat_p(lchnk, icol)
      loncol = get_lon_p(lchnk, icol)

      group = .false.
      mtapes = get_mtapes()
      ntapes = 0
      
!     check to be sure sites are listed on only one tape
      do it=1,mtapes
         if (ngroup(it) .gt. 0) ntapes=ntapes+1
      enddo
      if(ntapes>1)then
         write(6,*)'sites should only be listed on one auxillary tape'
	 call endrun()
      endif

!     look for grid cell in list of selected columns

      do it=1,mtapes
         do g=1,ngroup(it)
	    lon = tape(it)%column(g)%columnlon(1)
	    lat = tape(it)%column(g)%columnlat(1)
	    if(lon.ne.tape(it)%column(g)%columnlon(2).or.
     &	       lat.ne.tape(it)%column(g)%columnlat(2))then
	       write(6,*)'CRM column history error.'
	       write(6,*) 
     &	       'Cant handle more than one grid cell per column group'
	       write(6,*)
     &	       'Because dont know which column is last in group'
	       write(6,*)
     &	       'Also, columns in group could be on different nodes'
	       write(6,*)'lon,lat=',
     &             (tape(it)%column(g)%columnlon(i),i=1,2),
     &             (tape(it)%column(g)%columnlat(i),i=1,2)
               call endrun()
	    endif
	    if(lon.eq.loncol.and.lat.eq.latcol)then
	       ! have found a column to write
	       group = .true.
	       group_lon = lon
	       group_lat = lat
	       group_tape = it
	       group_number = g
            endif
         end do
      enddo
      
      if(.not.group)return ! grid cell is not in list of columns
      
!     grid cell is in list of columns
      
      it = group_tape
      g = group_number
      
!      write(6,*)'column found: lchnk,icol,group_lon,group_lat=',
!     &	   lchnk,icol,group_lon,group_lat,
!     &     ' group_tape,group_number=',it,g,' nfid=',nfid(g,it)

      if(nfid(g,it).eq.0)then

!        create new history file

         string =trim(interpret_filename_spec( hfilename_spec(it),
     &            (it-1),.false. ) )
         lenstring=len_trim(string)
         do i=1,lenstring-2
            if ( string(i:i+2) == ".nc" )then
	       nhfil(g,it) = string(1:i-1) // "." //
     &         trim(tape(it)%column(g)%lon_name) //
     &         "_" // trim(tape(it)%column(g)%lat_name)// ".nc"
            endif
	 enddo

         write(6,*)'Creating CRM netcdf history file ',
     &         trim(nhfil(g,it))
         call wrap_create (nhfil(g,it),nf_clobber, nfid(g,it))
         write(6,*)'nfid=',nfid(g,it)
         ret = nf_def_dim (nfid(g,it), 'time', nf_unlimited, timdim)
         ret = nf_def_dim (nfid(g,it), 'tbnd', 2, tbnddim)
         ret = nf_def_dim (nfid(g,it), 'chars', 8, chardim)
         ret = nf_def_dim (nfid(g,it), 'crm_x', crm_nx, crm_x_dim)
         ret = nf_def_dim (nfid(g,it), 'crm_y', crm_ny, crm_y_dim)
         ret = nf_def_dim (nfid(g,it), 'crm_z', crm_nz, crm_z_dim)
!	 write(6,*)tape(it)%column(g)%lat_name,tape(it)%column(g)%lon_name
         ret = nf_def_dim (nfid(g,it), 'lat',
     & 	       tape(it)%column(g)%num_lats, grouplatdim(g))
         call wrap_def_var (nfid(g,it),'lat',
     &	      nf_double,1,grouplatdim(g),glatvar(g))
         call wrap_put_att_text (nfid(g,it), glatvar(g),'long_name',
     &	   'latitude')
         call wrap_put_att_text (nfid(g,it), glatvar(g),'units',
     &	      'degrees_north')
         ret = nf_def_dim (nfid(g,it), 'lon',
     &	      tape(it)%column(g)%num_lons, grouplondim(g))
         call wrap_def_var (nfid(g,it),'lon',
     & 	      nf_double,1,grouplondim(g),glonvar(g))
         call wrap_put_att_text (nfid(g,it), glonvar(g),'long_name',
     &	      'longitude')
         call wrap_put_att_text (nfid(g,it), glonvar(g),'units',
     &	      'degrees_east')
         dimen1(1) = timdim
         dimen2t(1) = tbnddim
         dimen2t(2) = timdim

         dimen5g(1) = grouplondim(g)
         dimen5g(2) = grouplatdim(g)
         dimen5g(3) = crm_x_dim
         dimen5g(4) = crm_y_dim
         dimen5g(5) = timdim
	 dimen6g(1) = grouplondim(g)
	 dimen6g(2) = grouplatdim(g)
	 dimen6g(3) = crm_x_dim
	 dimen6g(4) = crm_y_dim
	 dimen6g(5) = crm_z_dim
	 dimen6g(6) = timdim
	 call wrap_def_var (nfid(g,it),'crm_x',nf_float,
     &	        1,crm_x_dim,crm_x_var)
         call wrap_put_att_text (nfid(g,it), crm_x_var, 'long_name',
     &	        'crm horiz. distance in x')
         call wrap_put_att_text (nfid(g,it), crm_x_var, 'units', 'km')

         call wrap_def_var (nfid(g,it),'crm_y',nf_float,
     &	        1,crm_y_dim,crm_y_var)
         call wrap_put_att_text (nfid(g,it), crm_y_var, 'long_name',
     &	         'crm horiz. distance in y')
         call wrap_put_att_text (nfid(g,it), crm_y_var, 'units', 'km')

         call wrap_def_var (nfid(g,it),'height',nf_float,
     &	        1, crm_z_dim, z_var)
         call wrap_put_att_text (nfid(g,it), z_var, 'long_name',
     &	        'height above surface')
         call wrap_put_att_text (nfid(g,it), z_var, 'units', 'm')

         call wrap_def_var (nfid(g,it),'crm_z',nf_float,
     &	        1, crm_z_dim, pres_var)
         call wrap_put_att_text (nfid(g,it), pres_var, 'long_name',
     &	        'crm nominal pressure')
         call wrap_put_att_text (nfid(g,it), pres_var, 'units', 'hPa')

         call wrap_def_var (nfid(g,it),'zs',nf_double, 0,0, zs_var)
         call wrap_put_att_text (nfid(g,it), zs_var, 'long_name',
     &	                        'surface elevation')
         call wrap_put_att_text (nfid(g,it), zs_var, 'units', 'm')

         call wrap_def_var (nfid(g,it),'ps',nf_double, 1,timdim,
     &                      ps_var(g,it))
         call wrap_put_att_text (nfid(g,it), ps_var(g,it), 'long_name',
     &	                        'surface pressure')
         call wrap_put_att_text (nfid(g,it), ps_var(g,it), 'units', 'hPa')

         call wrap_def_var (nfid(g,it),'ps0',nf_double, 0,0,ps0_var)
         call wrap_put_att_text (nfid(g,it), ps0_var, 'long_name',
     &	                        'reference pressure')
         call wrap_put_att_text (nfid(g,it), ps0_var, 'units', 'hPa')

         call get_ref_date(iyr, imon, iday, nbsec)
         nbdate = iyr*10000 + imon*100 + iday
         call wrap_def_var (nfid(g,it),'time',nf_double,1,timdim,
     &	        timeid(g,it))
         call wrap_put_att_text (nfid(g,it), timeid(g,it), 'long_name',
     &	        'time')
         str = 'days since ' // date2yyyymmdd(nbdate) // ' ' //
     &	        sec2hms(nbsec)
         call wrap_put_att_text (nfid(g,it), timeid(g,it), 'units', str)
         call wrap_put_att_text (nfid(g,it), timeid(g,it), 'calendar',
     &	        'noleap')

         call wrap_def_var (nfid(g,it),'time_bnds',nf_double,2,dimen2t,
     &	        tbndid(g,it))
         call wrap_put_att_text (nfid(g,it), tbndid(g,it), 'long_name',
     &	        'time interval endpoints')
         str = 'days since ' // date2yyyymmdd(nbdate) // ' ' //
     &	        sec2hms(nbsec)
         call wrap_put_att_text (nfid(g,it), tbndid(g,it), 'units', str)

!
! Character
!
         dimenchar(1) = chardim
         dimenchar(2) = timdim
         call wrap_def_var (nfid(g,it),'date_written',NF_CHAR,2,
     &	        dimenchar, date_writtenid(g,it))
         call wrap_def_var (nfid(g,it),'time_written',NF_CHAR,2,
     &	        dimenchar, time_writtenid(g,it))
         call wrap_def_var (nfid(g,it),'ndbase',NF_INT,0,0,
     &  	 ndbaseid(g,it))
         str = 'base day'
         call wrap_put_att_text (nfid(g,it), ndbaseid(g,it),
     &	        'long_name',str)

         call wrap_def_var (nfid(g,it),'nsbase',NF_INT,0,0,
     &          nsbaseid(g,it))
         str = 'seconds of base day'
         call wrap_put_att_text (nfid(g,it), nsbaseid(g,it),
     &	        'long_name',str)

         call wrap_def_var (nfid(g,it),'nbdate',NF_INT,0,0,
     &	      nbdateid(g,it))
         str = 'base date (YYYYMMDD)'
         call wrap_put_att_text (nfid(g,it), nbdateid(g,it),
     &	        'long_name',str)

         call wrap_def_var (nfid(g,it),'nbsec',NF_INT,0,0,nbsecid(g,it))
         str = 'seconds of base date'
         call wrap_put_att_text (nfid(g,it), nbsecid(g,it), 'long_name',
     &	        str)

         call wrap_def_var (nfid(g,it),'mdt',NF_INT,0,0,mdtid)
         call wrap_put_att_text (nfid(g,it), mdtid, 'long_name',
     &	        'timestep')
         call wrap_put_att_text (nfid(g,it), mdtid, 'units', 's')

!
! Character header information
!
         str = 'CF-1.0'
         call wrap_put_att_text (nfid(g,it), NF_GLOBAL, 'Conventions',
     &	        str)
         call wrap_put_att_text (nfid(g,it), NF_GLOBAL, 'source',
     &	        'CAM-CRM')
         logname = ' '
         call getenv ('LOGNAME',logname)
         if (logname=='        ') then
            write(6,*)'CRM: Cannot find LOGNAME environment variable'
            call endrun
         end if
         host = ' '
         call getenv ('HOST',host)
         call wrap_put_att_text (nfid(g,it), NF_GLOBAL, 'case',caseid)
         call wrap_put_att_text (nfid(g,it), NF_GLOBAL, 'title',ctitle)
         call wrap_put_att_text (nfid(g,it), NF_GLOBAL, 'logname',
     &	        logname)
         call wrap_put_att_text (nfid(g,it), NF_GLOBAL, 'host', host)
         call wrap_put_att_text (nfid(g,it), NF_GLOBAL, 'Version',
     &           '$Name:  $')
!
! Create variables for model timing and header information
!
         call wrap_def_var (nfid(g,it),'ndcur   ',nf_int,1,dimen1,
     &	        ndcurid(g,it))
         str = 'current day (from base day)'
         call wrap_put_att_text (nfid(g,it), ndcurid(g,it), 'long_name',
     &	        str)

         call wrap_def_var (nfid(g,it),'nscur   ',nf_int,1,dimen1,
     &	        nscurid(g,it))
         str = 'current seconds of current day'
         call wrap_put_att_text (nfid(g,it), nscurid(g,it), 'long_name',
     &	        str)

         call wrap_def_var (nfid(g,it),'date    ',nf_int,1,dimen1,
     &	        dateid(g,it))
         str = 'current date (YYYYMMDD)'
         call wrap_put_att_text (nfid(g,it), dateid(g,it), 'long_name',
     &	        str)

         call wrap_def_var (nfid(g,it),'datesec ',nf_int,1,dimen1,
     &	        datesecid(g,it))
         str = 'current seconds of current date'
         call wrap_put_att_text (nfid(g,it), datesecid(g,it), 'long_name',
     &	        str)

!        attributes of predicted fields

         ncreal = nf_float

         call wrap_def_var(nfid(g,it), 'CRM_U', ncreal,6,dimen6g,
     &	        varid_u(g,it))
         str = 'm/s'
         call wrap_put_att_text (nfid(g,it), varid_u(g,it), 'units',
     &	        str)
         str = 'CRM u-velocity'
         call wrap_put_att_text (nfid(g,it), varid_u(g,it),
     & 	        'long_name', str)

         call wrap_def_var(nfid(g,it), 'CRM_V', ncreal,6,dimen6g,
     &	        varid_v(g,it))
         str = 'm/s'
         call wrap_put_att_text (nfid(g,it), varid_v(g,it), 'units',
     &	        str)
         str = 'CRM v-velocity'
         call wrap_put_att_text (nfid(g,it), varid_v(g,it),
     & 	        'long_name', str)

         call wrap_def_var(nfid(g,it), 'CRM_W', ncreal,6,dimen6g,
     &	        varid_w(g,it))
         str = 'm/s'
         call wrap_put_att_text (nfid(g,it), varid_w(g,it), 'units',
     &	        str)
         str = 'CRM w-velocity'
         call wrap_put_att_text (nfid(g,it), varid_w(g,it),
     & 	        'long_name', str)

         call wrap_def_var(nfid(g,it), 'CRM_TABS', ncreal,6,dimen6g,
     &	        varid_t(g,it))
         str = 'deg K'
         call wrap_put_att_text (nfid(g,it), varid_t(g,it), 'units',
     &	        str)
         str = 'CRM temperature'
         call wrap_put_att_text (nfid(g,it), varid_t(g,it),
     & 	        'long_name', str)

         call wrap_def_var(nfid(g,it), 'CRM_QV', ncreal,6,dimen6g,
     &	        varid_qv(g,it))
         str = 'g/kg'
         call wrap_put_att_text (nfid(g,it), varid_qv(g,it), 'units',
     &	        str)
         str = 'CRM water vapor mixing ratio'
         call wrap_put_att_text (nfid(g,it), varid_qv(g,it),
     & 	        'long_name', str)

         call wrap_def_var(nfid(g,it), 'CRM_QC', ncreal,6,dimen6g,
     &	        varid_qc(g,it))
         str = 'g/kg'
         call wrap_put_att_text (nfid(g,it), varid_qc(g,it), 'units',
     &	        str)
         str = 'CRM cloud water mixing ratio'
         call wrap_put_att_text (nfid(g,it), varid_qc(g,it),
     & 	        'long_name', str)

         call wrap_def_var(nfid(g,it), 'CRM_QI', ncreal,6,dimen6g,
     &	        varid_qi(g,it))
         str = 'g/kg'
         call wrap_put_att_text (nfid(g,it), varid_qi(g,it), 'units',
     &	        str)
         str = 'CRM cloud ice mixing ratio'
         call wrap_put_att_text (nfid(g,it), varid_qi(g,it),
     & 	        'long_name', str)

         call wrap_def_var(nfid(g,it), 'CRM_QPC', ncreal,6,dimen6g,
     &	        varid_qr(g,it))
         str = 'g/kg'
         call wrap_put_att_text (nfid(g,it), varid_qr(g,it), 'units',
     &	        str)
         str = 'CRM precipitating water mixing ratio'
         call wrap_put_att_text (nfid(g,it), varid_qr(g,it),
     & 	        'long_name', str)

         call wrap_def_var(nfid(g,it), 'CRM_QPI', ncreal,6,dimen6g,
     &	        varid_qs(g,it))
         str = 'g/kg'
         call wrap_put_att_text (nfid(g,it), varid_qs(g,it), 'units',
     &	        str)
         str = 'CRM precipitating ice (snow+graupel) mixing ratio'
         call wrap_put_att_text (nfid(g,it), varid_qs(g,it),
     & 	        'long_name', str)

         call wrap_def_var(nfid(g,it), 'CRM_QRS', ncreal,6,dimen6g,
     &	        varid_qrs(g,it))
         str = 'Deg/s'
         call wrap_put_att_text (nfid(g,it), varid_qrs(g,it), 'units',
     &	        str)
         str = 'CRM sw radiative heating rate'
         call wrap_put_att_text (nfid(g,it), varid_qrs(g,it),
     & 	        'long_name', str)

         call wrap_def_var(nfid(g,it), 'CRM_QRL', ncreal,6,dimen6g,
     &	        varid_qrl(g,it))
         str = 'Deg/s'
         call wrap_put_att_text (nfid(g,it), varid_qrl(g,it), 'units',
     &	        str)
         str = 'CRM lw radiative heating rate'
         call wrap_put_att_text (nfid(g,it), varid_qrl(g,it),
     & 	        'long_name', str)

         call wrap_def_var(nfid(g,it), 'PRECSFC', ncreal,5,dimen5g,
     &          varid_prec(g,it))
         str = 'm/s'
         call wrap_put_att_text (nfid(g,it), varid_prec(g,it), 'units',
     &	        str)
         str = 'CRM precipitation'
         call wrap_put_att_text (nfid(g,it), varid_prec(g,it),
     & 	        'long_name', str)

         call wrap_def_var(nfid(g,it), 'CRM_FSDS', ncreal,5,dimen5g,
     &          varid_fsds(g,it))
         str = 'W/m2'
         call wrap_put_att_text (nfid(g,it), varid_fsds(g,it), 'units',
     &	        str)
         str = 'CRM Downward Solar at Surface'
         call wrap_put_att_text (nfid(g,it), varid_fsds(g,it),
     & 	        'long_name', str)

         call wrap_def_var(nfid(g,it), 'CRM_FSNS', ncreal,5,dimen5g,
     &          varid_fsns(g,it))
         str = 'W/m2'
         call wrap_put_att_text (nfid(g,it), varid_fsns(g,it), 'units',
     &	        str)
         str = 'CRM Net Downward Solar at Surface'
         call wrap_put_att_text (nfid(g,it), varid_fsns(g,it),
     & 	        'long_name', str)

         call wrap_def_var(nfid(g,it), 'CRM_FSDSC', ncreal,5,dimen5g,
     &          varid_fsdsc(g,it))
         str = 'W/m2'
         call wrap_put_att_text (nfid(g,it), varid_fsdsc(g,it), 'units',
     &	        str)
         str = 'CRM Clearsky Downward Solar at Surface'
         call wrap_put_att_text (nfid(g,it), varid_fsdsc(g,it),
     & 	        'long_name', str)

         call wrap_def_var(nfid(g,it), 'CRM_FSNT', ncreal,5,dimen5g,
     &          varid_fsnt(g,it))
         str = 'W/m2'
         call wrap_put_att_text (nfid(g,it), varid_fsnt(g,it), 'units',
     &	        str)
         str = 'CRM Net Upward Solar at Top of Atmosphere'
         call wrap_put_att_text (nfid(g,it), varid_fsnt(g,it),
     & 	        'long_name', str)

         call wrap_def_var(nfid(g,it), 'CRM_FSNTC', ncreal,5,dimen5g,
     &          varid_fsntc(g,it))
         str = 'W/m2'
         call wrap_put_att_text (nfid(g,it), varid_fsntc(g,it), 'units',
     &	        str)
         str = 'CRM Clearsky Net Upward Solar at Top of Atmosphere'
         call wrap_put_att_text (nfid(g,it), varid_fsntc(g,it),
     & 	        'long_name', str)

         call wrap_def_var(nfid(g,it), 'CRM_FSUT', ncreal,5,dimen5g,
     &          varid_fsut(g,it))
         str = 'W/m2'
         call wrap_put_att_text (nfid(g,it), varid_fsut(g,it), 'units',
     &	        str)
         str = 'CRM Upward Solar at Top of Atmosphere'
         call wrap_put_att_text (nfid(g,it), varid_fsut(g,it),
     & 	        'long_name', str)

         call wrap_def_var(nfid(g,it), 'CRM_FLDS', ncreal,5,dimen5g,
     &          varid_flds(g,it))
         str = 'W/m2'
         call wrap_put_att_text (nfid(g,it), varid_flds(g,it), 'units',
     &	        str)
         str = 'CRM Downward Longwave at Surface'
         call wrap_put_att_text (nfid(g,it), varid_flds(g,it),
     & 	        'long_name', str)

         call wrap_def_var(nfid(g,it), 'CRM_FLNS', ncreal,5,dimen5g,
     &          varid_flns(g,it))
         str = 'W/m2'
         call wrap_put_att_text (nfid(g,it), varid_flns(g,it), 'units',
     &	        str)
         str = 'CRM Net Upward Longwave at Surface'
         call wrap_put_att_text (nfid(g,it), varid_flns(g,it),
     & 	        'long_name', str)

         call wrap_def_var(nfid(g,it), 'CRM_FLNSC', ncreal,5,dimen5g,
     &          varid_flnsc(g,it))
         str = 'W/m2'
         call wrap_put_att_text (nfid(g,it), varid_flnsc(g,it), 'units',
     &	        str)
         str = 'CRM Clearsky Net Upward Longwave at Surface'
         call wrap_put_att_text (nfid(g,it), varid_flnsc(g,it),
     & 	        'long_name', str)

         call wrap_def_var(nfid(g,it), 'CRM_FLUT', ncreal,5,dimen5g,
     &          varid_flut(g,it))
         str = 'W/m2'
         call wrap_put_att_text (nfid(g,it), varid_flut(g,it), 'units',
     &	        str)
         str = 'CRM Upward Longwave at Top of Model'
         call wrap_put_att_text (nfid(g,it), varid_flut(g,it),
     & 	        'long_name', str)

         call wrap_def_var(nfid(g,it), 'CRM_FLUTC', ncreal,5,dimen5g,
     &          varid_flutc(g,it))
         str = 'W/m2'
         call wrap_put_att_text (nfid(g,it), varid_flutc(g,it), 'units',
     &	        str)
         str = 'CRM Clearsky Upward Longwave at Top of Model'
         call wrap_put_att_text (nfid(g,it), varid_flutc(g,it),
     & 	        'long_name', str)


!        end define mode

         ret = nf_enddef(nfid(g,it))
         write(6,*)'Successfully created CRM netcdf file ',nhfil(g,it)

         do i=1,crm_nx
            crm_x(i) = (i-1)*crm_dx*0.001
         end do
         call wrap_put_var_realx (nfid(g,it), crm_x_var, crm_x)
         do i=1,crm_ny
            crm_y(i) = (i-1)*crm_dy*0.001
         end do
         do i=1,crm_nz
           presref(i) = 0.01*ps0*(hyam(plev-i+1) + hybm(plev-i+1))
         end do
         call wrap_put_var_realx (nfid(g,it), crm_y_var, crm_y)
         call wrap_put_var_realx (nfid(g,it), z_var, z)
         call wrap_put_var_realx (nfid(g,it), pres_var, presref)
         call wrap_put_var_realx (nfid(g,it), zs_var, zs)
         call wrap_put_var_realx (nfid(g,it), ps0_var, ps0)
         call wrap_put_var_int (nfid(g,it), mdtid, dt)
         call wrap_put_var_realx (nfid(g,it), glonvar(g),
     &	      londeg(tape(it)%column(g)%columnlon(1):
     &               tape(it)%column(g)%columnlon(2),1))
         call wrap_put_var_realx (nfid(g,it), glatvar(g),
     &	      latdeg(tape(it)%column(g)%columnlat(1):
     &               tape(it)%column(g)%columnlat(2)))
!
! Model date info
!
         call wrap_put_var_int (nfid(g,it), ndbaseid(g,it), ndbase)
         call wrap_put_var_int (nfid(g,it), nsbaseid(g,it), nsbase)

         call wrap_put_var_int (nfid(g,it), nbdateid(g,it), nbdate)
         call wrap_put_var_int (nfid(g,it), nbsecid(g,it), nbsec)

	 call get_curr_time(ndcur, nscur)
         nscur=nscur-get_step_size() ! current time at beginning of timestep
         beg_time(g,it) = ndcur + nscur/86400._r8

         nfils(g,it) = 0
      
      endif

!     now write history data

      nfils(g,it) = nfils(g,it) + 1

!     call get_curr_date(iyr, imon, iday, ncsec)
      call get_prev_date(iyr, imon, iday, ncsec) ! current time at beginning of timestep
      ncdate = iyr*10000 + imon*100 + iday
      call get_curr_time(ndcur, nscur)
      nscur=nscur-get_step_size() ! current time at beginning of timestep

      nscur=nscur+time_crm ! add crm time to CAM time
      if(nscur<0)then
         nscur=nscur+86400
	 ndcur=ndcur-1
      endif
      ncsec=ncsec+time_crm ! add crm time to CAM time
!     write(6,*)'writing history at time ',ndcur,nscur
!     write(6,*)'for file ',nhfil(g,it)

      start(1)=nfils(g,it)
      count(1)=1
      call wrap_put_vara_int (nfid(g,it), ndcurid(g,it),start, count,
     &      ndcur)
      call wrap_put_vara_int (nfid(g,it), nscurid(g,it),start, count,
     &      nscur)
      call wrap_put_vara_int (nfid(g,it), dateid(g,it),start, count,
     &      ncdate)
      call wrap_put_vara_int (nfid(g,it), datesecid(g,it),start,count,
     &      ncsec)
      htime = ndcur + nscur/86400._r8
      call wrap_put_vara_realx (nfid(g,it), timeid(g,it), start, count,
     &      htime)
      call wrap_put_vara_realx (nfid(g,it), ps_var(g,it), start, count, ps)

      start(1) = 1
      start(2) = nfils(g,it)
      count(1) = 2
      count(2) = 1
      tdata(1) = beg_time(g,it)
      tdata(2) = htime
      call wrap_put_vara_realx (nfid(g,it), tbndid(g,it), start, count, tdata)
      beg_time(g,it) = htime  ! update beginning time of next interval

      start(1) = 1
      start(2) = nfils(g,it)
      count(1) = 8
      count(2) = 1
      call datetime (cdate, ctime)
      call wrap_put_vara_text (nfid(g,it), date_writtenid(g,it), start, 
     &       count, cdate)
      call wrap_put_vara_text (nfid(g,it), time_writtenid(g,it), start, 
     &       count, ctime)



      start(1) = group_lon-tape(it)%column(g)%columnlon(1)+1
      start(2) = group_lat-tape(it)%column(g)%columnlat(1)+1
      start(3) = 1
      start(4) = 1
      count(1) = 1
      count(2) = 1
      count(3) = crm_nx
      count(4) = crm_ny

!     3-D fields

      start(5) = 1
      start(6) = nfils(g,it)
      count(5) = crm_nz
      count(6) = 1

      do k=1,crm_nz
      do j=1,crm_ny
      do i=1,crm_nx
         work3d(i,j,k)=u(i,j,k)
      enddo
      enddo
      enddo
      call wrap_put_vara_real (nfid(g,it), varid_u(g,it),
     &      start, count, work3d)

      do k=1,crm_nz
      do j=1,crm_ny
      do i=1,crm_nx
         work3d(i,j,k)=v(i,j,k)
      enddo
      enddo
      enddo
      call wrap_put_vara_real (nfid(g,it), varid_v(g,it),
     &      start, count, work3d)

      do k=1,crm_nz
      do j=1,crm_ny
      do i=1,crm_nx
         work3d(i,j,k)=w(i,j,k)
      enddo
      enddo
      enddo
      call wrap_put_vara_real (nfid(g,it), varid_w(g,it),
     &      start, count, work3d)

      do k=1,crm_nz
      do j=1,crm_ny
      do i=1,crm_nx
      work3d(i,j,k)=tabs(i,j,k)
      enddo
      enddo
      enddo
      call wrap_put_vara_real (nfid(g,it), varid_t(g,it),
     &      start, count, work3d)

      do k=1,crm_nz
      do j=1,crm_ny
      do i=1,crm_nx
         work3d(i,j,k)=1000.*(q(i,j,k)-qn(i,j,k))
      enddo
      enddo
      enddo
      call wrap_put_vara_real (nfid(g,it), varid_qv(g,it),
     &      start, count, work3d)

	do k=1,crm_nz
	do j=1,crm_ny
	do i=1,crm_nx
	   work3d(i,j,k)=1000.*qn(i,j,k)*omegan(i,j,k)
	enddo
	enddo
	enddo
        call wrap_put_vara_real (nfid(g,it), varid_qc(g,it),
     &      start, count, work3d)

	do k=1,crm_nz
	do j=1,crm_ny
	do i=1,crm_nx
	   work3d(i,j,k)=1000.*qn(i,j,k)*(1.-omegan(i,j,k))
	enddo
	enddo
	enddo
        call wrap_put_vara_real (nfid(g,it), varid_qi(g,it),
     &      start, count, work3d)

	do k=1,crm_nz
	do j=1,crm_ny
	do i=1,crm_nx
	   work3d(i,j,k)=1000.*qp(i,j,k)*omegap(i,j,k)
	enddo
	enddo
	enddo
        call wrap_put_vara_real (nfid(g,it), varid_qr(g,it),
     &      start, count, work3d)

	do k=1,crm_nz
	do j=1,crm_ny
	do i=1,crm_nx
	   work3d(i,j,k)=1000.*qp(i,j,k)*(1.-omegap(i,j,k))
	enddo
	enddo
	enddo
        call wrap_put_vara_real (nfid(g,it), varid_qs(g,it),
     &      start, count, work3d)

	do k=1,crm_nz
	do j=1,crm_ny
	do i=1,crm_nx
	   work3d(i,j,k)=qrs_crm(i,j,k)
	enddo
	enddo
	enddo
        call wrap_put_vara_real (nfid(g,it), varid_qrs(g,it),
     &      start, count, work3d)

	do k=1,crm_nz
	do j=1,crm_ny
	do i=1,crm_nx
	   work3d(i,j,k)=qrl_crm(i,j,k)
	enddo
	enddo
	enddo
        call wrap_put_vara_real (nfid(g,it), varid_qrl(g,it),
     &      start, count, work3d)

!       2-D  fields

        start(5) = nfils(g,it)
        count(5) = 1

	do j=1,crm_ny
	do i=1,crm_nx
	   work2d(i,j)=prechist(i,j)
	enddo
	enddo
        call wrap_put_vara_real (nfid(g,it), varid_prec(g,it),
     &	    start, count, work2d)

	do j=1,crm_ny
	do i=1,crm_nx
	   work2d(i,j)=fsds(i,j)
	enddo
	enddo
        call wrap_put_vara_real (nfid(g,it), varid_fsds(g,it),
     &	    start, count, work2d)

	do j=1,crm_ny
	do i=1,crm_nx
	   work2d(i,j)=fsns(i,j)
	enddo
	enddo
        call wrap_put_vara_real (nfid(g,it), varid_fsns(g,it),
     &	    start, count, work2d)

	do j=1,crm_ny
	do i=1,crm_nx
	   work2d(i,j)=fsdsc(i,j)
	enddo
	enddo
        call wrap_put_vara_real (nfid(g,it), varid_fsdsc(g,it),
     &	    start, count, work2d)

	do j=1,crm_ny
	do i=1,crm_nx
	   work2d(i,j)=fsnt(i,j)
	enddo
	enddo
        call wrap_put_vara_real (nfid(g,it), varid_fsnt(g,it),
     &	    start, count, work2d)

	do j=1,crm_ny
	do i=1,crm_nx
	   work2d(i,j)=fsntc(i,j)
	enddo
	enddo
        call wrap_put_vara_real (nfid(g,it), varid_fsntc(g,it),
     &	    start, count, work2d)

	do j=1,crm_ny
	do i=1,crm_nx
	   work2d(i,j)=fsut(i,j)
	enddo
	enddo
        call wrap_put_vara_real (nfid(g,it), varid_fsut(g,it),
     &	    start, count, work2d)

	do j=1,crm_ny
	do i=1,crm_nx
	   work2d(i,j)=flds(i,j)
	enddo
	enddo
        call wrap_put_vara_real (nfid(g,it), varid_flds(g,it),
     &	    start, count, work2d)

	do j=1,crm_ny
	do i=1,crm_nx
	   work2d(i,j)=flns(i,j)
	enddo
	enddo
        call wrap_put_vara_real (nfid(g,it), varid_flns(g,it),
     &	    start, count, work2d)

	do j=1,crm_ny
	do i=1,crm_nx
	   work2d(i,j)=flnsc(i,j)
	enddo
	enddo
        call wrap_put_vara_real (nfid(g,it), varid_flnsc(g,it),
     &	    start, count, work2d)

	do j=1,crm_ny
	do i=1,crm_nx
	   work2d(i,j)=flut(i,j)
	enddo
	enddo
        call wrap_put_vara_real (nfid(g,it), varid_flut(g,it),
     &	    start, count, work2d)

	do j=1,crm_ny
	do i=1,crm_nx
	   work2d(i,j)=flutc(i,j)
	enddo
	enddo
        call wrap_put_vara_real (nfid(g,it), varid_flutc(g,it),
     &	    start, count, work2d)

!       now close if it's time

        if (nstep.ne.nstop) return !check to see if done with CRM loop
	if(ncsec.ge.86400.or.is_last_step())then ! check for end of day or job
!          close file on new day

           write(6,*)'closing CRM netcdf history file ',nhfil(g,it)
           ret = nf_close (nfid(g,it))
	   nfid(g,it)=0
        endif

	return
	end
