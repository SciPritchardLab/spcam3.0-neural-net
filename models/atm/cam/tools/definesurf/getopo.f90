subroutine getopo(ntopo, ntlon, ntlat, ftopo, htopo )

  use shr_kind_mod, only: r8 => shr_kind_r8

!
!       Read in /DSS/K0240K Navy 10 minute elevation dataset.
!       See DSS website info on DSS754.0 for more information.
!
!-----------------------------------------------------------------------
  implicit none
!-----------------------------------------------------------------------

  integer ntopo              ! Input:  Unit number of file
  integer ntlon              ! Input:  Number of lons for topo data
  integer ntlat              ! Input:  Number of lats for topo data

  real(r8)  ftopo(ntlon,ntlat) ! Output: Topo 0/1 flag array
  real(r8)  htopo(ntlon,ntlat) ! Output: Topographic heights (m)

!--------------------------Local variables------------------------------

  integer i,j,ii,jj,jt          ! Indices
  integer ilnd,iocn             ! counters for surrounding points
  integer lon,lat,in(9,30,30)   ! Input navy 10' by 10' data array
  real(r8) ft2m
  data ft2m /.3048/
!
! Dynamic
!
  integer mask(ntlon,ntlat)     ! 2-D land characteristics array
!-----------------------------------------------------------------------
! Read in topodata from S to N
! Input topographic height is in 100s of feet, convert to m upon input

  print *,' getopo: read TOPO data. ',ntlon,ntlat
  do j=1,36
    print *,j
    do i=1,72
      read (ntopo,'(8x,2i6)') lat,lon
      read(ntopo,1003)in
      do jj=1,30
        jt=(j-1)*30+jj
        do ii=1,30
          htopo((i-1)*30+ii,jt)=in(3,ii,jj)*100.*ft2m
          mask((i-1)*30+ii,jt) =in(6,ii,jj)
        enddo
      enddo
    enddo
  enddo
1003 format(5(1x,2i2,3i3,2i2,2i3))

!  Quality control: The documenation for this dataset defines values of
!  0-9 as different land types and 62 as ocean.  There are numerous
!  points with values which are greater than 9 and not 62.  
!  If one if these undocumented points is found, check the values of all
!  valid neighboring points and if half or more of them are land, make
!  the undocumented point land.
!  If all surrounding points are bad, make it land and print location.

  do j=1,ntlat
    do i=1,ntlon
      if (mask(i,j).gt.9 .and. mask(i,j).ne.62) then
        ilnd=0
        iocn=0
        do jj=max(1,j-1),min(ntlat,j+1)
          do ii=max(1,i-1),min(ntlon,i+1)
            if (mask(ii,jj).lt.9 .and. mask(ii,jj).ge.0) then
              ilnd=ilnd+1
            else if (mask(ii,jj).eq.62) then
              iocn=iocn+1
            endif
          enddo
        enddo
        if (ilnd .ge. iocn) then
          if (ilnd.eq.0) then
            print *,'rdtopo1: Bad data all around ',i,j,'  Using: ',mask(i,j-1)
            mask(i,j)=mask(i,j-1)
          else
            mask(i,j) = 1   ! For efficiency, make it land type 1
          endif
        else
          mask(i,j)=62
        endif
      endif
    enddo
  enddo
!
!  Set ftopo = 0 for ocean and 1 for land
!
  where (mask(:,:).eq.62)
    ftopo(:,:) = 0.0
  elsewhere
    ftopo(:,:) = 1.0
  end where
!
! Make the Caspian sea appear as ocean points
!     
  do j=1,ntlat
    do i=1,ntlon
      if ( (i.ge.280 .and. i.le.325)  .and. &
           (j.ge.760 .and. j.le.822)  .and. &
           (mask(i,j).eq.4         )        )   ftopo(i,j)=0.0
    enddo
  enddo
  print *,' getopo:  TOPO data read.'
  return
end subroutine getopo

