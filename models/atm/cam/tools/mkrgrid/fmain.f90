program fmain
  use constants
  use control
!
! $Id: fmain.f90,v 1.2 2000/05/18 19:35:41 rosinski Exp $
!
  implicit none
!
! Local workspace
!
  character*80 arg, file1, file2
  character*256 :: cmdline             ! input command line
  character*8 cvsid
  data cvsid/'none'/

  integer n
  integer nargs
      
  integer iargc
  external iargc
!
! Default settings before parsing argument list
!
  file1 = ' '
  file2 = ' '
  default_interp = 'linear'
  monosiz = 0
  verbose = .false.
  silent = .false.
  reverse = .false.
  iceflag = .false.

  nargs = iargc()
  n = 1
  cmdline = char(10) // 'mkrgrid '
  do while (n .le. nargs)
    arg = ' '
    call getarg (n, arg)
    n = n + 1
    select case (arg)
    case ('-c')
      default_interp = 'cubic'
      cmdline = trim(cmdline) // ' -c'
    case ('-f')
      default_interp = 'fourier'
      cmdline = trim(cmdline) // ' -f'
    case ('-i')
      iceflag = .true.
      cmdline = trim(cmdline) // ' -i'
    case ('-l')
      default_interp = 'linear'
      cmdline = trim(cmdline) // ' -l'
    case ('-m')
      call getarg (n, arg)
      n = n + 1
      monosiz = monosiz + 1
      monolist(monosiz) = arg(1:8)
      cmdline = trim(cmdline) // ' -m ' // trim(monolist(monosiz))
    case ('-r')
      reverse = .true.
      cmdline = trim(cmdline) // ' -r'
    case ('-s')
      silent = .true.
      cmdline = trim(cmdline) // ' -s'
    case ('-v')
      verbose = .true.
      cmdline = trim(cmdline) // ' -v'
    case default
      if (file1 .eq. ' ') then
        file1 = arg
      else if (file2 .eq. ' ') then
        file2 = arg
      else
        write (6,*) 'Argument ', arg,' is not known'
        call usage_exit (' ')
      end if
      cmdline = trim(cmdline) // ' ' // trim(arg)
    end select
  end do
  
  if (file1.eq.' ' .or. file2.eq.' ') then
    call usage_exit ('Must enter an input file and an output file')
  else if (silent .and. verbose) then
    call usage_exit ('-s cannot be specified with -v')
  else if (iceflag .and. reverse) then
    call usage_exit ('Cannot convert ice flags Reduced -> full ')
  end if

  call mkrgrid (cvsid, file1, file2, cmdline)
                
  stop
end program

subroutine usage_exit (arg)
  implicit none
  character*(*) arg

  if (arg.ne.' ') write (6,*) arg
  write (6,*) 'Usage: mkrgrid [-c] [-f] [-i] [-l] [-m field]... ', &
              '[-r] [-s] [-v] infile outfile [ < namelist]'
  stop 999
end subroutine
