module constants
    use precision_kinds, only: prec

        implicit none
    
    private

    real(prec), parameter, public :: pi = acos(-1.0_prec)

end module constants
