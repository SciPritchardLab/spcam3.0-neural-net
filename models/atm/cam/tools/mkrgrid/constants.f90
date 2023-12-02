module constants
  use shr_kind_mod, only: r8 => shr_kind_r8

  real(r8), parameter :: fillvalue  = 1.d36

  integer, parameter  :: maxnlat    = 10000
  integer, parameter  :: maxmonosiz = 1000
end module
