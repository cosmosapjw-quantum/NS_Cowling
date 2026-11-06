! Floating-point representation only; physical dimensions are not Fortran types.
module ns_precision
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    private
    integer, parameter, public :: wp = real64
end module ns_precision
