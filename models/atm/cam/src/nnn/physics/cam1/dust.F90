#include <misc.h>
#include <params.h>

module dust

implicit none

  private          ! Make default type private to the module
  save
!
! Public interfaces
!
  public dust_number,ndst,sz_nbr,dst_src_nbr                        ! number of dust constituents

  integer, parameter:: nx_dust = 4
  integer, parameter:: ndst =4
  integer, parameter:: dst_src_nbr =3
  integer, parameter:: sz_nbr =200
contains

  function dust_number()
    implicit none
    integer dust_number
    dust_number = nx_dust
  end function dust_number

end module dust
