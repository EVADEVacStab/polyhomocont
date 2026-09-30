module test_system__closed_form
!! A user-supplied system in closed form, as a user of the library would
!! write it:
!!   F_1 = x1^2 + x2^2 - 5
!!   F_2 = x1*x2^2 - 2*x2 + x1
!! with degrees (2, 3), together with its homogenization written out by
!! hand for comparison with the default `evaluate_homogeneous`.

  use polyhomocont, only : wp
  use polyhomocont, only : poly_system

  implicit none

  private

  public :: closed_form_system
  public :: homogenized_by_hand

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

  subroutine homogenized_by_hand(z, f, jac)
    !! F^h_1 = z1^2 + z2^2 - 5 w^2, F^h_2 = z1 z2^2 - 2 z2 w^2 + z1 w^2.

    complex(wp), intent(in) :: z(3)
    complex(wp), intent(out) :: f(2)
    complex(wp), intent(out) :: jac(2, 3)

    f(1) = z(1)**2 + z(2)**2 - 5.0_wp*z(3)**2
    f(2) = z(1)*z(2)**2 - 2.0_wp*z(2)*z(3)**2 + z(1)*z(3)**2
    jac(1, :) = [2.0_wp*z(1), 2.0_wp*z(2), -10.0_wp*z(3)]
    jac(2, :) = [z(2)**2 + z(3)**2,  &
      2.0_wp*z(1)*z(2) - 2.0_wp*z(3)**2,  &
      -4.0_wp*z(2)*z(3) + 2.0_wp*z(1)*z(3)]

  end subroutine homogenized_by_hand

end module test_system__closed_form

program test_system
!! Tests of the abstract `poly_system` with a user-supplied closed-form
!! system: Jacobian against finite differences and the default
!! homogenization against the hand-written homogenized system.

  use polyhomocont, only : wp
  use polyhomocont__linalg, only : norm_inf
  use test_system__closed_form, only : closed_form_system
  use test_system__closed_form, only : homogenized_by_hand
  use testing, only : assert_small
  use testing, only : assert_equal
  use testing, only : jacobian_fd_error
  use testing, only : homogeneous_jacobian_fd_error
  use testing, only : homogeneity_error
  use testing, only : finish

  implicit none

  type(closed_form_system) :: sys
  complex(wp) :: x(2)
  complex(wp) :: z(3)
  complex(wp) :: f(2)
  complex(wp) :: fref(2)
  complex(wp) :: jac(2, 3)
  complex(wp) :: jref(2, 3)
  integer, allocatable :: d(:)

  call assert_equal("nvars", sys%nvars(), 2)
  d = sys%degrees()
  call assert_equal("degree 1", d(1), 2)
  call assert_equal("degree 2", d(2), 3)

  x = [(0.3_wp, -1.2_wp), (2.1_wp, 0.4_wp)]
  call assert_small("affine Jacobian vs finite differences",  &
    jacobian_fd_error(sys, x), 1.0e-8_wp)

  ! Default homogenization at generic points, including a point with
  ! small w (close to infinity)
  z = [(0.3_wp, -1.2_wp), (2.1_wp, 0.4_wp), (0.7_wp, 0.2_wp)]
  call sys%evaluate_homogeneous(z, f, jac)
  call homogenized_by_hand(z, fref, jref)
  call assert_small("F^h vs hand-written", norm_inf(f - fref), 1.0e-13_wp)
  call assert_small("dF^h vs hand-written",  &
    maxval(abs(jac - jref)), 1.0e-13_wp)
  call assert_small("dF^h vs finite differences",  &
    homogeneous_jacobian_fd_error(sys, z), 1.0e-8_wp)
  call assert_small("homogeneity",  &
    homogeneity_error(sys, z, (0.6_wp, -1.3_wp)), 1.0e-13_wp)

  ! Close to infinity, dF^h/dw suffers from cancellation in
  ! d_i F_i - sum_j x_j dF_i/dx_j (relative error ~ eps/|w|)
  z(3) = (1.0e-6_wp, -2.0e-6_wp)
  call sys%evaluate_homogeneous(z, f, jac)
  call homogenized_by_hand(z, fref, jref)
  call assert_small("F^h vs hand-written, small w",  &
    norm_inf(f - fref)/norm_inf(fref), 1.0e-9_wp)
  call assert_small("dF^h vs hand-written, small w",  &
    maxval(abs(jac - jref))/maxval(abs(jref)), 1.0e-8_wp)

  call finish()

end program test_system
