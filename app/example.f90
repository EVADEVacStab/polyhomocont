module example__user_system
!! A system supplied through its own evaluator instead of the sparse
!! representation: the intersection of the circle x^2 + y^2 = 4 with the
!! cubic y = x^3 - 2x + 1/2.

  use polyhomocont, only : wp
  use polyhomocont, only : poly_system

  implicit none

  private

  public :: circle_cubic

  type, extends(poly_system) :: circle_cubic
  contains
    procedure :: nvars => cc_nvars
    procedure :: degrees => cc_degrees
    procedure :: evaluate => cc_evaluate
  end type circle_cubic

contains

  function cc_nvars(self) result(n)

    class(circle_cubic), intent(in) :: self
    integer :: n

    n = 2

  end function cc_nvars

  function cc_degrees(self) result(d)

    class(circle_cubic), intent(in) :: self
    integer, allocatable :: d(:)

    d = [2, 3]

  end function cc_degrees

  subroutine cc_evaluate(self, x, f, jac)

    class(circle_cubic), intent(inout) :: self
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(out) :: f(:)
    complex(wp), intent(out) :: jac(:, :)

    f(1) = x(1)**2 + x(2)**2 - 4.0_wp
    f(2) = x(1)**3 - 2.0_wp*x(1) + 0.5_wp - x(2)
    jac(1, :) = [2.0_wp*x(1), 2.0_wp*x(2)]
    jac(2, :) = [3.0_wp*x(1)**2 - 2.0_wp, (-1.0_wp, 0.0_wp)]

  end subroutine cc_evaluate

end module example__user_system

program example
!! Solves two small systems with polyhomocont: one given in the sparse
!! representation, one through a user-supplied evaluator.

  use polyhomocont, only : wp
  use polyhomocont, only : status_ok
  use polyhomocont, only : sparse_system
  use polyhomocont, only : solve
  use polyhomocont, only : solve_result
  use example__user_system, only : circle_cubic

  implicit none

  type(sparse_system) :: sys
  type(circle_cubic) :: user
  type(solve_result) :: res
  integer :: status

  ! x^2 + y^2 = 5, x y = 2. Each column of the exponent matrix is the
  ! exponent vector of one term.
  call sys%init(2)
  call sys%add_equation([1.0_wp, 1.0_wp, -5.0_wp],  &
    reshape([2, 0,  0, 2,  0, 0], [2, 3]), status)
  call sys%add_equation([1.0_wp, -2.0_wp],  &
    reshape([1, 1,  0, 0], [2, 2]), status)
  call solve(sys, res, status)
  if (status /= status_ok) error stop "solve failed"
  print "(a)", "x^2 + y^2 = 5, x y = 2:"
  call report(res)

  call solve(user, res, status)
  if (status /= status_ok) error stop "solve failed"
  print "(/,a)", "x^2 + y^2 = 4, y = x^3 - 2x + 1/2 (user evaluator):"
  call report(res)

contains

  subroutine report(res)

    type(solve_result), intent(in) :: res

    integer :: j

    print "(a,i0,a,i0,a,i0,a,i0,a,i0)", "  paths: ", res%npaths,  &
      ", solutions: ", res%nsolutions, " (real: ", res%nreal,  &
      "), at infinity: ", res%ninfinite, ", failed: ", res%nfailed
    do j = 1, res%nsolutions
      print "(2x,a,2(f10.6,sp,f10.6,ss,a))", "(",  &
        res%solutions(1, j), "i, ", res%solutions(2, j), "i)"
    end do
    print "(a)", "  real solutions:"
    do j = 1, res%nreal
      print "(2x,2f12.8)", res%real_solutions(:, j)
    end do

  end subroutine report

end program example
