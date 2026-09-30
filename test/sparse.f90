program test_sparse
!! Tests of the sparse polynomial system: values and Jacobians against
!! hand-written expressions and finite differences, exact homogenization
!! (also at infinity), random sparse systems and input validation.

  use polyhomocont, only : wp
  use polyhomocont, only : status_ok
  use polyhomocont, only : err_invalid_input
  use polyhomocont, only : sparse_system
  use polyhomocont__random, only : rng
  use polyhomocont__linalg, only : norm_inf
  use testing, only : assert_status
  use testing, only : assert_small
  use testing, only : assert_equal
  use testing, only : assert_true
  use testing, only : jacobian_fd_error
  use testing, only : homogeneous_jacobian_fd_error
  use testing, only : homogeneity_error
  use testing, only : finish

  implicit none

  type(rng) :: gen

  call gen%seed(11)
  call check_small_system()
  call check_random_systems()
  call check_validation()
  call finish()

contains

  subroutine check_small_system()
    !! F_1 = x1^2 + x2^2 - 5, F_2 = x1*x2^2 - 2*x2 + x1 (the closed-form
    !! system of test/system.f90), with an extra zero term in F_2.

    type(sparse_system) :: sys
    complex(wp) :: x(2)
    complex(wp) :: z(3)
    complex(wp) :: f(2)
    complex(wp) :: fref(2)
    complex(wp) :: jac(2, 2)
    complex(wp) :: jref(2, 2)
    integer, allocatable :: d(:)
    integer :: status

    call sys%init(2)
    call sys%add_equation([1.0_wp, 1.0_wp, -5.0_wp],  &
      reshape([2, 0,  0, 2,  0, 0], [2, 3]), status)
    call assert_status("small: add equation 1", status)
    call sys%add_equation(  &
      [(1.0_wp, 0.0_wp), (-2.0_wp, 0.0_wp), (0.0_wp, 0.0_wp),  &
       (1.0_wp, 0.0_wp)],  &
      reshape([1, 2,  0, 1,  5, 5,  1, 0], [2, 4]), status)
    call assert_status("small: add equation 2", status)

    call assert_equal("small: nvars", sys%nvars(), 2)
    call assert_equal("small: nequations", sys%nequations(), 2)
    call assert_equal("small: zero term dropped", sys%nterms(2), 3)
    d = sys%degrees()
    call assert_true("small: degrees (2, 3)", all(d == [2, 3]))

    x = [(0.3_wp, -1.2_wp), (2.1_wp, 0.4_wp)]
    call sys%evaluate(x, f, jac)
    fref(1) = x(1)**2 + x(2)**2 - 5.0_wp
    fref(2) = x(1)*x(2)**2 - 2.0_wp*x(2) + x(1)
    jref(1, :) = [2.0_wp*x(1), 2.0_wp*x(2)]
    jref(2, :) = [x(2)**2 + 1.0_wp, 2.0_wp*x(1)*x(2) - 2.0_wp]
    call assert_small("small: F", norm_inf(f - fref), 1.0e-14_wp)
    call assert_small("small: dF/dx", maxval(abs(jac - jref)), 1.0e-14_wp)

    ! Homogenized system at a finite point and exactly at infinity
    z = [(0.3_wp, -1.2_wp), (2.1_wp, 0.4_wp), (0.7_wp, 0.2_wp)]
    call check_homogenized(sys, z, "small: finite point")
    z(3) = (0.0_wp, 0.0_wp)
    call check_homogenized(sys, z, "small: at infinity")

    ! Points with zero coordinates
    x = [(0.0_wp, 0.0_wp), (1.5_wp, -0.5_wp)]
    call assert_small("small: Jacobian at x1 = 0",  &
      jacobian_fd_error(sys, x), 1.0e-8_wp)
    x = [(0.0_wp, 0.0_wp), (0.0_wp, 0.0_wp)]
    call sys%evaluate(x, f, jac)
    jref(1, :) = [(0.0_wp, 0.0_wp), (0.0_wp, 0.0_wp)]
    jref(2, :) = [(1.0_wp, 0.0_wp), (-2.0_wp, 0.0_wp)]
    call assert_small("small: Jacobian at x = 0",  &
      maxval(abs(jac - jref)), 0.0_wp)

  end subroutine check_small_system

  subroutine check_homogenized(sys, z, name)
    !! Compares the homogenized small system with the hand-written
    !! F^h_1 = z1^2 + z2^2 - 5 w^2, F^h_2 = z1 z2^2 - 2 z2 w^2 + z1 w^2.

    type(sparse_system), intent(inout) :: sys
    complex(wp), intent(in) :: z(3)
    character(len=*), intent(in) :: name

    complex(wp) :: f(2)
    complex(wp) :: fref(2)
    complex(wp) :: jh(2, 3)
    complex(wp) :: jhref(2, 3)

    call sys%evaluate_homogeneous(z, f, jh)
    fref(1) = z(1)**2 + z(2)**2 - 5.0_wp*z(3)**2
    fref(2) = z(1)*z(2)**2 - 2.0_wp*z(2)*z(3)**2 + z(1)*z(3)**2
    jhref(1, :) = [2.0_wp*z(1), 2.0_wp*z(2), -10.0_wp*z(3)]
    jhref(2, :) = [z(2)**2 + z(3)**2,  &
      2.0_wp*z(1)*z(2) - 2.0_wp*z(3)**2,  &
      -4.0_wp*z(2)*z(3) + 2.0_wp*z(1)*z(3)]
    call assert_small(name//": F^h", norm_inf(f - fref), 1.0e-14_wp)
    call assert_small(name//": dF^h", maxval(abs(jh - jhref)),  &
      1.0e-14_wp)

  end subroutine check_homogenized

  subroutine check_random_systems()
    !! Random sparse systems in 4 variables with random exponents (0..3)
    !! and complex coefficients: Jacobians against finite differences and
    !! homogeneity of the homogenized system.

    integer, parameter :: n = 4
    integer, parameter :: nt = 6
    type(sparse_system) :: sys
    complex(wp) :: coefs(nt)
    integer :: exps(n, nt)
    complex(wp) :: x(n)
    complex(wp) :: z(n + 1)
    integer :: status
    integer :: trial
    integer :: i
    integer :: j
    integer :: k
    real(wp) :: err_j
    real(wp) :: err_jh
    real(wp) :: err_h

    err_j = 0.0_wp
    err_jh = 0.0_wp
    err_h = 0.0_wp
    do trial = 1, 5
      call sys%init(n)
      do i = 1, n
        do k = 1, nt
          coefs(k) = gen%normal_complex()
          do j = 1, n
            exps(j, k) = int(4.0_wp*gen%uniform())
          end do
        end do
        call sys%add_equation(coefs, exps, status)
        call assert_status("random: add equation", status)
      end do
      do j = 1, n
        x(j) = gen%normal_complex()
        z(j) = gen%normal_complex()
      end do
      z(n + 1) = gen%normal_complex()
      err_j = max(err_j, jacobian_fd_error(sys, x))
      err_jh = max(err_jh, homogeneous_jacobian_fd_error(sys, z))
      err_h = max(err_h, homogeneity_error(sys, z, gen%normal_complex()))
    end do
    call assert_small("random: Jacobian vs finite differences", err_j,  &
      1.0e-7_wp)
    call assert_small("random: homogenized Jacobian vs FD", err_jh,  &
      1.0e-7_wp)
    call assert_small("random: homogeneity", err_h, 1.0e-12_wp)

  end subroutine check_random_systems

  subroutine check_validation()
    !! Invalid input must be rejected with err_invalid_input.

    type(sparse_system) :: sys
    integer :: status

    call sys%init(2)
    call sys%add_equation([1.0_wp, 2.0_wp],  &
      reshape([1, 0, 0], [3, 1]), status)
    call assert_true("validation: wrong exponent rows",  &
      status == err_invalid_input)
    call sys%add_equation([1.0_wp, 2.0_wp],  &
      reshape([1, 0], [2, 1]), status)
    call assert_true("validation: wrong number of terms",  &
      status == err_invalid_input)
    call sys%add_equation([1.0_wp],  &
      reshape([1, -1], [2, 1]), status)
    call assert_true("validation: negative exponent",  &
      status == err_invalid_input)
    call sys%add_equation([3.0_wp, 0.0_wp],  &
      reshape([0, 0,  1, 1], [2, 2]), status)
    call assert_true("validation: constant equation",  &
      status == err_invalid_input)
    call assert_equal("validation: nothing added", sys%nequations(), 0)

    call sys%add_equation([1.0_wp], reshape([1, 0], [2, 1]), status)
    call assert_true("validation: first equation ok", status == status_ok)
    call sys%add_equation([1.0_wp], reshape([0, 1], [2, 1]), status)
    call assert_true("validation: second equation ok",  &
      status == status_ok)
    call sys%add_equation([1.0_wp], reshape([1, 1], [2, 1]), status)
    call assert_true("validation: too many equations",  &
      status == err_invalid_input)

  end subroutine check_validation

end program test_sparse
