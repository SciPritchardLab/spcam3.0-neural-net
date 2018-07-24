module nt_FunctionsModule
    public :: nt_logistic, logisticd, const, constd

    contains
    
        !Logistic activation function a/(b+exp(-x))
        !args - array of parameter:
        !args(1) - a
        !args(2) - b

        function nt_logistic(x, args) result(fx)
            real, intent(in) :: args(0:)
            real, intent(in) :: x
            real :: fx

            fx = args(1) / (args(2) + exp(-x))

        end function nt_logistic

        !Logistic function derivative
        function logisticd(args, x) result(fx)
            implicit none
            real, intent(in) :: args(:)
            real, intent(in) :: x
            real :: fx,expmx,d

            expmx = exp(-x)
            d     = args(2)+expmx
            fx    = args(1)*expmx/d/d

        end function logisticd

        function const(args, x) result(fx)
            implicit none
            real, intent(in) :: args(:)
            real, intent(in) :: x
            real :: fx

            fx=args(1)

        end function const

        function constd(args, x) result(fx)
            implicit none
            real, intent(in) :: args(:)
            real, intent(in) :: x
            real :: fx

            fx=0

        end function constd


        function nt_logsig(x) result(fx)
            real, intent(in) :: x
            real :: fx

            fx = 1. / (1. + exp(-x))

        end function nt_logsig

        function logsigd(x) result(fx)
            implicit none
            real, intent(in) :: x
            real :: fx,expmx,d

            expmx = exp(-x)
            d     = 1.+expmx
            fx    = expmx/d/d

        end function logsigd

        function nt_tansig(x) result(fx)
            real, intent(in) :: x
            real :: fx

            fx = 2./(1.+exp(-2.*n))-1.

        end function nt_tansig

        function tansigd(x) result(fx)
            implicit none
            real, intent(in) :: x
            real :: fx,expmx,d

            expmx = exp(-2.*x)
            d     = 1.+expmx
            fx    = 4.*expmx/d/d

        end function tansigd




end module nt_FunctionsModule
