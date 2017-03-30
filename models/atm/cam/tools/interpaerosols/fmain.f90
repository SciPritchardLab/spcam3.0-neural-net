program fmain
!------------------------------------------------------------------
! Purpose:
!   Usage: interpaerosols [-d] [-v] <infile.nc> <outfile.nc>
!
!   Horizontally interpolate file containing aerosol masses (infile.nc)
!   to CAM grid (as specified by outfile.nc).  Perl script REGRID.pl
!   contained in this directory can be used to create a template
!   outfile.nc.
!
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
   use globals

   implicit none

   include 'netcdf.inc'
!
! Local workspace
!
   character(len=80) :: arg = ' '        ! cmd line arg
   character(len=256) :: infile = ' '    ! input file name
   character(len=256) :: outfile = ' '   ! output file name
   character(len=256) :: cmdline = ' '   ! command line

   integer :: n                          ! argument counter
   integer :: nargs                      ! number of command line arguments
   integer :: ncprec = nf_float          ! default precision to write data

   logical :: verbose = .false.          ! verbose output

   integer iargc
   external iargc
!
! Default settings before parsing argument list
!
   nargs = iargc()
   n = 1
   cmdline = 'interpaerosols '
   do while (n <= nargs)
      arg = ' '
      call getarg (n, arg)
      n = n + 1
      
      select case (arg)
      case ('-d')
         ncprec = nf_double
         cmdline = trim(cmdline) // ' -d'
      case ('-v')
         verbose = .true.
         cmdline = trim(cmdline) // ' -v'
      case default
         if (infile == ' ') then
            infile = arg
         else if (outfile == ' ') then
            outfile = arg
         else
            write (6,*) 'Argument ', arg,' is not known'
            call usage_exit (' ')
         end if
         cmdline = trim(cmdline) // ' ' // trim(arg)
      end select
   end do
  
   if (infile == ' ' .or. outfile == ' ') then
      call usage_exit ('Must enter an input file and an output file')
   end if
!
! Open input and output netcdf files
!
   call wrap_open (infile, NF_NOWRITE, ncidi)
   call wrap_open (outfile, NF_WRITE, ncido)
   call wrap_redef (ncido)                      ! start in define mode on output file
!
! Get input and output dimensions for x, y, time
! Input file first
!
   call wrap_inq_dimid (ncidi, 'lon', londimidi)
   call wrap_inq_dimlen (ncidi, londimidi, nxi)

   call wrap_inq_dimid (ncidi, 'lat', latdimidi)
   call wrap_inq_dimlen (ncidi, latdimidi, nyi)
   
   call wrap_inq_dimid (ncidi, 'lev', levdimidi)
   call wrap_inq_dimlen (ncidi, levdimidi, nz)        ! nz is the same for input and output

   call wrap_inq_dimid (ncidi, 'time', timedimidi)
   call wrap_inq_dimlen (ncidi, timedimidi, ntime)
!
! Ensure ntime is 12.  Reason is because mean-preserving code for the moment assumes it.
!
   if (ntime /= 12) then
      write(6,*) 'Size of input time dimension must be 12'
      stop 999
   end if
!
! Now output file
!
   call wrap_inq_dimid (ncido, 'lon', londimido)
   call wrap_inq_dimlen (ncido, londimido, nxo)

   call wrap_inq_dimid (ncido, 'lat', latdimido)
   call wrap_inq_dimlen (ncido, latdimido, nyo)
!
! Assume z dimensions on output file don't exist (code will crash if they do)
! This way it will be guaranteed to match input file.
!
   call wrap_def_dim (ncido, 'lev', nz, levdimido)
   call wrap_def_dim (ncido, 'ilev', nz+1, ilevdimido)
!
! If a time dimension is not found on the output file, create it
!
   if (nf_inq_dimid (ncido, 'time', timedimido) /= nf_noerr) then
      call wrap_def_dim (ncido, 'time', nf_unlimited, timedimido)
   end if
!
! Add global attributes
!
   call addglobal (ncido, cmdline)
!
! Call driver code to read the input data and do the interpolations and/or copies
!
   call driver (ncprec)
!
! Close netcdf files.  This is crucial for output file in particular in order
! to ensure that all data get written to the file.
!
   call wrap_close (ncidi)
   call wrap_close (ncido)

   write(6,*)'Successfully interpolated aerosol data to file ', trim (outfile)

   stop 0
end program fmain

subroutine usage_exit (arg)
   implicit none
   character*(*) arg

   if (arg /= ' ') write (6,*) arg
   write (6,*) 'Usage: interpaerosols [-d] [-v] infile.nc outfile.nc'
   stop 999
end subroutine usage_exit
