module polyhomocont__config
!! Working precision, constants and status codes shared by all
!! polyhomocont modules.

  use, intrinsic :: iso_fortran_env, only : real64

  implicit none

  private

  integer, parameter, public :: wp = real64
    !! Working precision.

  real(wp), parameter, public :: pi = 3.14159265358979323846_wp

  ! Status codes of library procedures. `status_ok` (0) means success.
  integer, parameter, public :: status_ok = 0
  integer, parameter, public :: err_invalid_input = 1
    !! Inconsistent or invalid input (dimensions, degrees, ...).
  integer, parameter, public :: err_singular_matrix = 2
    !! A linear system could not be solved (exactly singular matrix).
  integer, parameter, public :: err_numerical = 3
    !! Generic numerical failure.

end module polyhomocont__config
