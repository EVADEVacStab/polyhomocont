module polyhomocont
!! Top-level module of the polyhomocont package: re-exports the public
!! API of the individual modules.
!!
!! polyhomocont computes all isolated complex solutions of square
!! systems of polynomial equations with homotopy continuation.

  use polyhomocont__config, only : wp
  use polyhomocont__config, only : status_ok
  use polyhomocont__config, only : err_invalid_input
  use polyhomocont__config, only : err_singular_matrix
  use polyhomocont__config, only : err_numerical
  use polyhomocont__system, only : poly_system
  use polyhomocont__sparse, only : sparse_system

  implicit none

  private

  ! Kind and status codes
  public :: wp
  public :: status_ok
  public :: err_invalid_input
  public :: err_singular_matrix
  public :: err_numerical
  ! Polynomial systems
  public :: poly_system
  public :: sparse_system

end module polyhomocont
