module precision_kinds
    use, intrinsic :: iso_fortran_env, only: real32, real64, real128
    
        implicit none
    
    private

    integer, parameter, public :: sp = real32
    integer, parameter, public :: dp = real64
    integer, parameter, public :: qp = real128
    
    integer, parameter, public :: prec = dp

end module precision_kinds
