module test_solve__closed_form
!! User-supplied closed-form version of the system
!!   F_1 = x1^2 + x2^2 - 5,  F_2 = x1*x2^2 - 2*x2 + x1,
!! to check that user evaluators give the same solutions as the sparse
!! representation.

  use polyhomocont, only : wp
  use polyhomocont, only : poly_system

  implicit none

  private

  public :: closed_form_system

  type, extends(poly_system) :: closed_form_system
  contains
    procedure :: nvars => cf_nvars
    procedure :: degrees => cf_degrees
    procedure :: evaluate => cf_evaluate
  end type closed_form_system

contains

  function cf_nvars(self) result(n)

    class(closed_form_system), intent(in) :: self
    integer :: n

    n = 2

  end function cf_nvars

  function cf_degrees(self) result(d)

    class(closed_form_system), intent(in) :: self
    integer, allocatable :: d(:)

    d = [2, 3]

  end function cf_degrees

  subroutine cf_evaluate(self, x, f, jac)

    class(closed_form_system), intent(inout) :: self
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(out) :: f(:)
    complex(wp), intent(out) :: jac(:, :)

    f(1) = x(1)**2 + x(2)**2 - 5.0_wp
    f(2) = x(1)*x(2)**2 - 2.0_wp*x(2) + x(1)
    jac(1, 1) = 2.0_wp*x(1)
    jac(1, 2) = 2.0_wp*x(2)
    jac(2, 1) = x(2)**2 + 1.0_wp
    jac(2, 2) = 2.0_wp*x(1)*x(2) - 2.0_wp

  end subroutine cf_evaluate

end module test_solve__closed_form

program test_solve
!! Tests of the solver on systems with known solutions: univariate
!! polynomials, a 2x2 system (independence of the random seed),
!! solutions at infinity, a double root, a user-supplied system, and
!! input validation.

  use polyhomocont, only : wp
  use polyhomocont, only : err_invalid_input
  use polyhomocont, only : sparse_system
  use polyhomocont, only : solve
  use polyhomocont, only : solve_options
  use polyhomocont, only : solve_result
  use polyhomocont__linalg, only : norm_inf
  use test_solve__closed_form, only : closed_form_system
  use testing, only : assert_status
  use testing, only : assert_small
  use testing, only : assert_equal
  use testing, only : assert_true
  use testing, only : set_distance
  use testing, only : finish

  implicit none

  call check_univariate()
  call check_2x2()
  call check_infinity()
  call check_double_root()
  call check_user_system()
  call check_invalid()
  call finish()

contains

  subroutine check_univariate()
    !! p(x) = (x - 1)(x + 2)(x - 3i)(x + 3i)(x - 1/2): 5 solutions, 3 real.

    complex(wp), parameter :: roots(5) = [(1.0_wp, 0.0_wp),  &
      (-2.0_wp, 0.0_wp), (0.0_wp, 3.0_wp), (0.0_wp, -3.0_wp),  &
      (0.5_wp, 0.0_wp)]
    type(sparse_system) :: sys
    type(solve_result) :: res
    complex(wp) :: c(0:5)
    integer :: status
    integer :: k

    ! Expand the product of (x - r_k): c(j) is the coefficient of x^j
    c = (0.0_wp, 0.0_wp)
    c(0) = (1.0_wp, 0.0_wp)
    do k = 1, 5
      c(1:k) = c(0:k-1) - roots(k)*c(1:k)
      c(0) = -roots(k)*c(0)
    end do

    call sys%init(1)
    call sys%add_equation(c, reshape([0, 1, 2, 3, 4, 5], [1, 6]), status)
    call solve(sys, res, status)
    call assert_status("univariate: solve", status)
    call assert_equal("univariate: paths", res%npaths, 5)
    call assert_equal("univariate: solutions", res%nsolutions, 5)
    call assert_equal("univariate: real solutions", res%nreal, 3)
    call assert_equal("univariate: failed", res%nfailed, 0)
    call assert_equal("univariate: duplicates", res%nduplicates, 0)
    call assert_small("univariate: distance to roots",  &
      set_distance(res%solutions, reshape(roots, [1, 5])), 1.0e-12_wp)

  end subroutine check_univariate

  subroutine check_2x2()
    !! x^2 + y^2 = 5, x y = 2: 4 real solutions, for several seeds.

    type(sparse_system) :: sys
    type(solve_result) :: res
    type(solve_options) :: opts
    complex(wp) :: sols(2, 4)
    real(wp) :: dist
    integer :: status
    integer :: seed
    integer :: nsol_ok
    integer :: nreal_ok

    sols = reshape([(1.0_wp, 0.0_wp), (2.0_wp, 0.0_wp),  &
      (2.0_wp, 0.0_wp), (1.0_wp, 0.0_wp),  &
      (-1.0_wp, 0.0_wp), (-2.0_wp, 0.0_wp),  &
      (-2.0_wp, 0.0_wp), (-1.0_wp, 0.0_wp)], [2, 4])
    call sys%init(2)
    call sys%add_equation([1.0_wp, 1.0_wp, -5.0_wp],  &
      reshape([2, 0,  0, 2,  0, 0], [2, 3]), status)
    call sys%add_equation([1.0_wp, -2.0_wp],  &
      reshape([1, 1,  0, 0], [2, 2]), status)

    dist = 0.0_wp
    nsol_ok = 0
    nreal_ok = 0
    do seed = 1, 10
      opts%seed = seed
      call solve(sys, res, status, opts)
      if (res%nsolutions == 4) nsol_ok = nsol_ok + 1
      if (res%nreal == 4) nreal_ok = nreal_ok + 1
      dist = max(dist, set_distance(res%solutions, sols))
    end do
    call assert_equal("2x2: 4 solutions for all 10 seeds", nsol_ok, 10)
    call assert_equal("2x2: 4 real solutions for all 10 seeds",  &
      nreal_ok, 10)
    call assert_small("2x2: distance to known solutions", dist,  &
      1.0e-12_wp)
    call assert_small("2x2: imaginary parts",  &
      maxval(abs(aimag(res%solutions))), 1.0e-12_wp)

  end subroutine check_2x2

  subroutine check_infinity()
    !! x y = 1, x^2 = 4: Bezout number 4, but only the 2 finite solutions
    !! (2, 1/2) and (-2, -1/2); the other paths go to the (singular)
    !! point (0 : 1 : 0) at infinity and must not be reported as
    !! solutions.

    type(sparse_system) :: sys
    type(solve_result) :: res
    complex(wp) :: sols(2, 2)
    integer :: status

    sols = reshape([(2.0_wp, 0.0_wp), (0.5_wp, 0.0_wp),  &
      (-2.0_wp, 0.0_wp), (-0.5_wp, 0.0_wp)], [2, 2])
    call sys%init(2)
    call sys%add_equation([1.0_wp, -1.0_wp],  &
      reshape([1, 1,  0, 0], [2, 2]), status)
    call sys%add_equation([1.0_wp, -4.0_wp],  &
      reshape([2, 0,  0, 0], [2, 2]), status)
    call solve(sys, res, status)
    call assert_status("infinity: solve", status)
    call assert_equal("infinity: paths", res%npaths, 4)
    call assert_equal("infinity: solutions", res%nsolutions, 2)
    call assert_small("infinity: distance to known solutions",  &
      set_distance(res%solutions, sols), 1.0e-12_wp)
    call assert_equal("infinity: remaining paths infinite/failed",  &
      res%ninfinite + res%nfailed + res%nsingular, 2)

  end subroutine check_infinity

  subroutine check_double_root()
    !! (x - 1)^2 (x + 1): the simple root -1 is the only nonsingular
    !! solution; the two paths to the double root 1 are singular or fail
    !! close to t = 0 (no end game yet).

    type(sparse_system) :: sys
    type(solve_result) :: res
    integer :: status

    call sys%init(1)
    call sys%add_equation([1.0_wp, -1.0_wp, -1.0_wp, 1.0_wp],  &
      reshape([3, 2, 1, 0], [1, 4]), status)
    call solve(sys, res, status)
    call assert_status("double root: solve", status)
    call assert_equal("double root: nonsingular solutions",  &
      res%nsolutions, 1)
    call assert_small("double root: solution -1",  &
      abs(res%solutions(1, 1) + 1.0_wp), 1.0e-12_wp)
    call assert_equal("double root: remaining paths singular/failed",  &
      res%nsingular + res%nfailed, 2)

  end subroutine check_double_root

  subroutine check_user_system()
    !! The closed-form user system and its sparse representation must
    !! give the same solutions, which must solve the system.

    type(closed_form_system) :: user
    type(sparse_system) :: sys
    type(solve_result) :: res_user
    type(solve_result) :: res_sparse
    complex(wp) :: f(2)
    complex(wp) :: jac(2, 2)
    real(wp) :: resid
    integer :: status
    integer :: j

    call sys%init(2)
    call sys%add_equation([1.0_wp, 1.0_wp, -5.0_wp],  &
      reshape([2, 0,  0, 2,  0, 0], [2, 3]), status)
    call sys%add_equation([1.0_wp, -2.0_wp, 1.0_wp],  &
      reshape([1, 2,  0, 1,  1, 0], [2, 3]), status)

    call solve(user, res_user, status)
    call assert_status("user system: solve", status)
    call solve(sys, res_sparse, status)
    call assert_status("user system: solve sparse", status)

    call assert_true("user system: solutions found",  &
      res_user%nsolutions > 0)
    call assert_equal("user system: same number of solutions",  &
      res_user%nsolutions, res_sparse%nsolutions)
    call assert_small("user system: same solutions",  &
      set_distance(res_user%solutions, res_sparse%solutions), 1.0e-11_wp)

    resid = 0.0_wp
    do j = 1, res_user%nsolutions
      call user%evaluate(res_user%solutions(:, j), f, jac)
      resid = max(resid, norm_inf(f)  &
        /max(1.0_wp, norm_inf(res_user%solutions(:, j)))**3)
    end do
    call assert_small("user system: residual", resid, 1.0e-13_wp)
    call assert_equal("user system: duplicates", res_user%nduplicates, 0)

  end subroutine check_user_system

  subroutine check_invalid()
    !! A system with fewer equations than variables is rejected.

    type(sparse_system) :: sys
    type(solve_result) :: res
    integer :: status

    call sys%init(2)
    call sys%add_equation([1.0_wp, -1.0_wp],  &
      reshape([2, 0,  0, 0], [2, 2]), status)
    call solve(sys, res, status)
    call assert_true("invalid: rejected", status == err_invalid_input)

  end subroutine check_invalid

end program test_solve
