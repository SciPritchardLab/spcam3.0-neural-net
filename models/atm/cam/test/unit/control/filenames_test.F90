module filenames_test
!
! Unit test for the filenames module. This uses the unit_driver
! program which expects two public interfaces:
!
! subroutine list_tests( tests, ntests )
! Where tests is a character array of ntests
!
! subroutine general_test
! General test.
!
! subroutine fail_test( test_name )
! Where test_name is one of the valid tests listed from filenames_list_tests.
!
  implicit none
!
! Public interfaces
!
  public list_tests    ! List tests that can be performed
  public general_test  ! General test
  public fail_test     ! Fail test
!
! Everything else private
!
  PRIVATE
  integer, parameter :: ntests = 9   ! number of tests
  ! Character description of tests
  character(len=180), parameter :: test_types(ntests) = (/ &
                         'general                                              ',  &
                         'failtest:init_filepaths:notabsdirname                ',  &
                         'failtest:get_archivedir:typewrong                    ',  &
                         'failtest:interpret_filename_spec:histtapenonotthere  ',  &
                         'failtest:interpret_filename_spec:nttoolarge          ',  &
                         'failtest:interpret_filename_spec:badexpandsion       ',  &
                         'failtest:interpret_filename_spec:filenametoolong     ',  &
                         'failtest:interpret_filename_spec:emptyhfilename_spec ',  &
                         'failtest:interpret_filename_spec:spaceinfilenamespec '   &
                                /)

CONTAINS

subroutine list_tests( tests, maxtest, ntest )
  character(len=*), intent(out) :: tests(*)   ! List of valid tests
  integer, intent(in) :: maxtest              ! Max number of tests
  integer, intent(out) :: ntest               ! Number of tests
!
! List tests that can be performed
!
  integer :: i

  ntest = ntests
  tests(:ntest) = test_types(:ntest)
  do i = 1, ntests
    write(6,'(a)') test_types(i)
  end do
end subroutine list_tests

subroutine init_hfilename_spec( ptapes, nhtfrq, hfilename_spec )
  integer, intent(in) ::  ptapes
  integer, intent(in) :: nhtfrq(ptapes)
  character(len=*), intent(out) ::  hfilename_spec(ptapes)

  integer nt
  do nt = 1, ptapes
     if ( len_trim(hfilename_spec(nt)) == 0 )then
        if ( nhtfrq(nt) == 0 )then
            hfilename_spec(nt) = '%c.cam2.h%t.%y-%m.nc'
        else
            hfilename_spec(nt) = '%c.cam2.h%t.%y-%m-%d-%s.nc'
        end if
     end if
  end do
end subroutine init_hfilename_spec

subroutine fail_test( test )
!
! Run a test that ensures a test that should fail works correctly
!
  use filenames, only: get_archivedir, &
                       interpret_filename_spec, &
                       init_filepaths
  character(len=*), intent(in) :: test
  integer, parameter :: nlen = 72
  character(len=300) :: rgpath
  character(len=nlen) :: bldfln
  integer  :: i

  print *, 'Fail test: ', test
  select case ( test )
     case( 'failtest:init_filepaths:notabsdirname' )
        call init_filepaths( archivedirname='subdir/dir' )
     case( 'failtest:get_archivedir:typewrong' )
        rgpath = get_archivedir( 'bad' )
     case( 'failtest:interpret_filename_spec:histtapenonotthere' )
        call init_filepaths(  )
        bldfln = interpret_filename_spec( 'history.nc' )
     case( 'failtest:interpret_filename_spec:nttoolarge' )
        call init_filepaths(  )
        bldfln = interpret_filename_spec( 'history%t.nc', number=10000 )
     case( 'failtest:interpret_filename_spec:badexpandsion' )
        call init_filepaths(  )
        bldfln = interpret_filename_spec( 'history%z.nc', number=1 )
     case( 'failtest:interpret_filename_spec:filenametoolong' )
        call init_filepaths(  )
        bldfln = interpret_filename_spec( &
              'history%y-%y-%y-%-%y-%y-%y-%s-%s-%s-%s-%s-%s-%s-%s-%s-%s-%s.nc', &
              number=1 )
     case( 'failtest:interpret_filename_spec:emptyhfilename_spec' )
        call init_filepaths(  )
        bldfln = interpret_filename_spec( '', number=1 )
     case( 'failtest:interpret_filename_spec:spaceinfilenamespec' )
        call init_filepaths( )
        bldfln = interpret_filename_spec( '%c %d.nc', number=0 )
     case default
        write(6,*) 'Test is not valid: ', test
  end select
end subroutine fail_test

subroutine general_test
!
! Put the filenames module through general testing
!
  use filenames, only: interpret_filename_spec, get_archivedir, &
                       init_filepaths, &
                       nrevsn, &
                       caseid, rest_pfile
  use shr_kind_mod, only: r8 => shr_kind_r8
  use time_manager, only: timemgr_init, advance_timestep, dtime, start_ymd, &
                          stop_ymd, is_last_step
#if ( defined SPMD )
   use mpishorthand, only: mpicom
#endif
! for nlend
#include <comctl.h>
  character(len=290) :: file, path
  character(len=290), parameter :: typename(3) = (/'init', 'rest', 'hist'/)
  character(len=300) rgpath       ! restart filepath
  integer nstep, nt
  integer, parameter :: ptapes = 6
  integer :: nhtfrq(ptapes)
  character(len=300) hfilename_spec(ptapes)
  character(len=300) adir
  integer, parameter :: p_spec_tests = 19
  character(len=300), parameter :: spec_tests(p_spec_tests) = (/    &
		'h%t%y%m%d.nc                                     ',  &
		'cam2.%t.%y.%m.%d.nc                              ',  &
		'cam2.%y.%m.%d.nc                                 ',  &
		'h.%s.%m.%d.%y.nc                                 ',  &
		'auxdir/h.%s.%m.%d.%y.nc                          ',  &
		'auxdir/%c/%y/h:%m.%d-%s.nc                       ',  &
		'ABShistory#$@!&*():;{}[]-_+=|<>,.?0123456789~`.nc',  &
		'Threepercentsigns:%%%%%%.nc                      ',  &
		'history.nc                                       ',  &
		'%c%y%m%d%s                                       ',  &
		'%c.%y.%m.%d.%s                                   ',  &
		'%c.:%y.:%m.:%d.:%s                               ',  &
		'%c.:-%y.:-%m.:-%d.:-%s                           ',  &
		'a%c.:-%y.:-%m.:-%d.:-%sa                         ',  &
		'ab%c.:-%y.:-%m.:-%d.:-%sab                       ',  &
		'abc%c.:-%y.:-%m.:-%d.:-%sabc                     ',  &
		'abcd%c.:-%y.:-%m.:-%d.:-%sabcd                   ',  &
		'abcde%c.:-%y.:-%m.:-%d.:-%sabcde                 ',  &
		'abcdef%c.:-%y.:-%m.:-%d.:-%sabcdef               '   &
		 /)
  ! test with: start_ymd = 150001019, nt=1
  character(len=72), parameter :: spec_answr(p_spec_tests) = (/    &
		'h0150001019.nc                                   ',  &
		'cam2.0.15000.10.19.nc                            ',  &
		'cam2.15000.10.19.nc                              ',  &
		'h.00000.10.19.15000.nc                           ',  &
		'auxdir/h.00000.10.19.15000.nc                    ',  &
		'auxdir/camTest/15000/h:10.19-00000.nc            ',  &
		'ABShistory#$@!&*():;{}[]-_+=|<>,.?0123456789~`.nc',  &
		'Threepercentsigns:%%%.nc                         ',  &
		'history.nc                                       ',  &
		'camTest15000101900000                            ',  &
		'camTest.15000.10.19.00000                        ',  &
		'camTest.:15000.:10.:19.:00000                    ',  &
		'camTest.:-15000.:-10.:-19.:-00000                ',  &
		'acamTest.:-15000.:-10.:-19.:-00000a              ',  &
		'abcamTest.:-15000.:-10.:-19.:-00000ab            ',  &
		'abccamTest.:-15000.:-10.:-19.:-00000abc          ',  &
		'abcdcamTest.:-15000.:-10.:-19.:-00000abcd        ',  &
		'abcdecamTest.:-15000.:-10.:-19.:-00000abcde      ',  &
		'abcdefcamTest.:-15000.:-10.:-19.:-00000abcdef    '   &
		 /)
  integer :: i   ! loop index
  integer, parameter :: ntape(6) = (/ 9, 10, 11, 99, 100, 101 /)
  logical :: first
!
! Filename specifiers
! (%c = caseid, $y = year, $m = month, $d = day, $s = seconds in day, %t = number)
!
   character(len=256) :: ifilename_spec = '%c.cam2.i.%y-%m-%d-%s'
   character(len=256) :: rfilename_spec = '%c.cam2.r.%y-%m-%d-%s'
   character(len=256) :: rafilename_spec = '%c.cam2.ra.%y-%m-%d-%s'
   character(len=256) :: rhfilename_spec = '%c.cam2.rh%t.%y-%m-%d-%s'


#if ( defined SPMD )
   call spmdinit ()
#endif
  caseid = "camTest"
  dtime = 3456    ! Use 3456 for a long time-step that gets irregular times
  start_ymd = 901
  stop_ymd = 250901
  nrefrq = 1
  call timemgr_init
  call init_filepaths( )
  write(6,*) 'Default: rest_pfile = ', rest_pfile
  hfilename_spec(:) = ' '
  nhtfrq(:) = 1
  nhtfrq(1) = 0
  call init_hfilename_spec( ptapes, nhtfrq, hfilename_spec )

  do i = 1, 3
    path = get_archivedir( typename(i) )
    print *, typename(i), ' path = ', path
    if ( trim(path) .ne. "/ERIK/csm/camTest/atm/"//trim(typename(i))//"/" )then
      write(6,*) 'ERROR:: bad path'
      call endrun
    end if
  end do
  rest_pfile = "/home/USER/rpointer/cam2.case.rpointer"
  adir = '/USER/csm/case'
  call init_filepaths( archivedirname=adir )
  print *, 'Test: rest_pfile = ', rest_pfile
  do i = 1, 3
    path = get_archivedir( typename(i) )
    print *, typename(i), ' path = ', path
    if ( trim(path) .ne. trim(adir)//"/"//trim(typename(i))//"/" )then
      write(6,*) 'ERROR:: bad path'
      call endrun
    end if
  end do
  file = interpret_filename_spec( rfilename_spec )
  print *, 'Test setting restart filepath with root directory: '
  rgpath = '/'//trim(file)
  print *, 'Now test setting stuff (as per a restart)'
  call init_filepaths( )
  do nt = 1, ptapes
    file = interpret_filename_spec( rhfilename_spec, number=(nt-1) )
    rgpath = get_archivedir('rest')//file
    print *, nt, ' history restart filepath = ', rgpath
  end do
  print *, 'Test that changing archive_dir is ok: '
  adir = '/CCSM/csm/case/cam2/data'
  call init_filepaths( archivedirname=adir );
  do i = 1, 3
    path = get_archivedir( typename(i) )
    print *, typename(i), ' path = ', path
    if ( trim(path) .ne. trim(adir)//"/"//trim(typename(i))//"/" )then
      write(6,*) 'ERROR:: bad path'
      call endrun
    end if
  end do
  print *, 'Test for a high year number and when primary files are not monthly:'
  start_ymd = 2147481230
  stop_ymd = start_ymd + 1
  print *, 'start_ymd = ', start_ymd
  call timemgr_init
  call init_filepaths( )
  nhtfrq(:) = 1
  call init_hfilename_spec( ptapes, nhtfrq, hfilename_spec )
  nt = 1
  file = interpret_filename_spec( hfilename_spec(nt), number=(nt-1) )
  write(6,'(a)') 'non-monthly primary history file: ', trim(file)
  print *, 'Another high year number test: '
  start_ymd = 150001019
  stop_ymd = start_ymd + 1
  print *, 'start_ymd = ', start_ymd
  call timemgr_init
  call init_filepaths( )
  file = interpret_filename_spec( hfilename_spec(nt), number=(nt-1) )
  print *, 'primary history file: ', trim(file)
  print *, 'Now test for different hfilename_spec: '
  do i = 1, p_spec_tests
    file = interpret_filename_spec( spec_tests(i), number=(nt-1) )
    print *, 'spec tests: ', trim(spec_tests(i)), ' file = ', file
    if ( trim(file) .ne. trim(spec_answr(i)) )then
      write(6,*) 'ERROR:: hfilename_spec test is wrong should = ', spec_answr(i)
      call endrun
    end if
  end do
!
! Test filenames over a longer time-series
!
  start_ymd = 99951019
  stop_ymd = 100050901
  call timemgr_init
  nlend = .false.
  first = .true.
  do while( .not. nlend )
    path = get_archivedir( 'rest' )
    print *, 'Build filenames return results: '
    file = interpret_filename_spec( ifilename_spec )
    print *, 'init file: ', trim(file)
    file = interpret_filename_spec( rfilename_spec )
    rgpath = trim(path) // trim(file)
    print *, 'master restart file: ', trim(file)
    if ( first )then
      file = interpret_filename_spec( rafilename_spec )
      print *, 'abs/ems restart file: ', trim(file)
      rgpath = trim(path) // trim(file)
      do nt = 1, ptapes
        print *, 'Tape number: ', nt
        file = interpret_filename_spec( rhfilename_spec, number=(nt-1) )
        rgpath = trim(path) // trim(file)
        print *, 'his restart file: ', trim(file)
        file = interpret_filename_spec( hfilename_spec(nt), number=(nt-1) )
        print *, 'history file: ', trim(file)
      end do
      do i = 1, size(ntape)
        file = interpret_filename_spec( hfilename_spec(2), number=(ntape(i)-1) )
        print *, 'history file: ', trim(file), ' nt = ', ntape(i)
      end do
    end if
    call advance_timestep
    if (is_last_step()) nlend = .true.
    first = .false.
  end do
#if ( defined SPMD )
   call mpibarrier (mpicom)
   call mpifinalize
#endif
end subroutine general_test

end module filenames_test

program unit_driver
!
! Driver for unit tests of the filenames module.
!
! Two interfaces are expected for unit-test modules:
!
! a list_tests, general_test, and fail_test subroutines.
!
! for example:
!
! subroutine list_tests( tests, ntests )
! Where tests is a character array of ntests
!
  use filenames_test, only: list_tests, general_test, fail_test
  implicit none
  character(len=72) :: test      ! Test to run
  integer ntest   ! Number of tests that can be done
  integer i       ! Loop index
  integer, parameter :: maxtests = 1000
  character(len=72) :: test_types(maxtests)

  print *, "filename_tests"
  print *, 'Enter a type of test to run: (list to list tests)'
  read(5,'(a)') test
  select case ( test )
     case( 'list' )
        call list_tests( test_types, maxtests, ntest )
     case( 'general' )
        call general_test
     case default
        call fail_test( test )
  end select
end program unit_driver
