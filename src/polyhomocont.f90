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
  use polyhomocont__tracker, only : tracker_options
  use polyhomocont__tracker, only : path_success
  use polyhomocont__tracker, only : path_failed_min_step
  use polyhomocont__tracker, only : path_failed_max_steps
  use polyhomocont__solve, only : solve
  use polyhomocont__solve, only : solve_options
  use polyhomocont__solve, only : solve_result
  use polyhomocont__solve, only : path_result

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
  ! Solver
  public :: solve
  public :: solve_options
  public :: solve_result
  public :: path_result
  public :: tracker_options
  public :: path_success
  public :: path_failed_min_step
  public :: path_failed_max_steps

end module polyhomocont
