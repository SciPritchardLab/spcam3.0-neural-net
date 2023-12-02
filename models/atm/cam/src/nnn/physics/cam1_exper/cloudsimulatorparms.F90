module cloudsimulatorparms
   use shr_kind_mod, only: r8 => shr_kind_r8
   
   integer, parameter :: ntau=7    ! number of optical depth ranges
   integer, parameter :: npres=7   ! number of pressure ranges

   logical :: doisccp = .false.    ! whether to do ISCCP calcs and I/O     

   real(r8) :: prlim(npres+1) = (/0., 180., 310., 440., 560., 680., 800., 1000./)
!
!JR The commented out line is what the documentation says
!JR The uncomment line is what the code does
!   real(r8) :: taulim(ntau+1)  = (/0., 0.1, 1.3, 3.6, 9.4, 23., 60., 379./)
   real(r8) :: taulim(ntau+1)  = (/0., 0.3, 1.3, 3.6, 9.4, 23., 60., 379./)
   
end module cloudsimulatorparms
