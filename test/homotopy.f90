program test_homotopy
!! Tests of the total-degree homotopy: start solutions, derivatives with
!! respect to z and t, and input validation.

  use polyhomocont, only : wp
  use polyhomocont, only : err_invalid_input
  use polyhomocont, only : sparse_system
  use polyhomocont__random, only : rng
  use polyhomocont__linalg, only : norm_inf
  use polyhomocont__homotopy, only : total_degree_homotopy
  use testing, only : assert_status
  use testing, only : assert_small
  use testing, only : assert_equal
  use testing, only : assert_true
  use testing, only : finish

  implicit none

  type(sparse_system), target :: sys
  type(sparse_system), target :: incomplete
  type(total_degree_homotopy) :: hom
  type(rng) :: gen
  integer :: status

  ! F_1 = x1^2 + x2^2 - 5, F_2 = x1*x2^2 - 2*x2 + x1; degrees (2, 3)
  call sys%init(2)
  call sys%add_equation([1.0_wp, 1.0_wp, -5.0_wp],  &
    reshape([2, 0,  0, 2,  0, 0], [2, 3]), status)
  call sys%add_equation([1.0_wp, -2.0_wp, 1.0_wp],  &
    reshape([1, 2,  0, 1,  1, 0], [2, 3]), status)

  call gen%seed(5)
  call hom%init(sys, gen, status)
  call assert_status("init", status)
  call assert_equal("nunknowns", hom%nunknowns(), 3)
  call assert_equal("nstart = Bezout number", hom%nstart, 6)
  call assert_small("|gamma| - 1", abs(hom%gamma) - 1.0_wp, 1.0e-15_wp)

  call check_start_solutions()

  ! Re-initialization (new gamma and patch) must work
  call hom%init(sys, gen, status)
  call assert_status("re-init", status)
  call check_start_solutions()
  call check_derivatives((0.37_wp, 0.0_wp))
  call check_derivatives((1.0_wp, 0.0_wp))
  call check_derivatives((0.0_wp, 0.0_wp))
  call check_derivatives((0.02_wp, -0.05_wp))

  ! A system with fewer equations than variables must be rejected
  call incomplete%init(2)
  call incomplete%add_equation([1.0_wp, -1.0_wp],  &
    reshape([2, 0,  0, 0], [2, 2]), status)
  call hom%init(incomplete, gen, status)
  call assert_true("init: non-square system rejected",  &
    status == err_invalid_input)

  call finish()

contains

  subroutine check_start_solutions()
    !! All start solutions solve H(z, 1) = 0 (including the patch) and are
    !! pairwise distinct.

    complex(wp) :: z(3, 6)
    complex(wp) :: h(3)
    complex(wp) :: hz(3, 3)
    complex(wp) :: ht(3)
    real(wp) :: res
    real(wp) :: dmin
    integer :: k
    integer :: l

    res = 0.0_wp
    do k = 1, hom%nstart
      call hom%start_solution(k, z(:, k))
      call hom%evaluate(z(:, k), (1.0_wp, 0.0_wp), h, hz, ht)
      ! Relative to the scale |z|^3 of the terms of the cubic equation
      res = max(res, norm_inf(h)/max(1.0_wp, norm_inf(z(:, k)))**3)
    end do
    call assert_small("start solutions: relative residual", res,  &
      1.0e-14_wp)

    dmin = huge(1.0_wp)
    do k = 1, hom%nstart
      do l = k + 1, hom%nstart
        dmin = min(dmin, norm_inf(z(:, k) - z(:, l)))
      end do
    end do
    call assert_true("start solutions: distinct", dmin > 1.0e-3_wp)

  end subroutine check_start_solutions

  subroutine check_derivatives(t)
    !! dH/dz and dH/dt against central finite differences at a random z
    !! (for complex t, the derivative along the real t direction).

    complex(wp), intent(in) :: t

    complex(wp) :: z(3)
    complex(wp) :: zp(3)
    complex(wp) :: h(3)
    complex(wp) :: hp(3)
    complex(wp) :: hm(3)
    complex(wp) :: hz(3, 3)
    complex(wp) :: ht(3)
    complex(wp) :: dz(3, 3)
    complex(wp) :: dt(3)
    real(wp), parameter :: step = 1.0e-5_wp
    real(wp) :: err
    character(len=48) :: name
    integer :: j

    do j = 1, 3
      z(j) = gen%normal_complex()
    end do
    call hom%evaluate(z, t, h, hz, ht)

    err = 0.0_wp
    do j = 1, 3
      zp = z
      zp(j) = z(j) + step
      call hom%evaluate(zp, t, hp, dz, dt)
      zp(j) = z(j) - step
      call hom%evaluate(zp, t, hm, dz, dt)
      err = max(err, norm_inf((hp - hm)/(2.0_wp*step) - hz(:, j)))
    end do
    write(name, "(a,f5.2,sp,f6.2,a)") "dH/dz vs FD, t = ", t, "i"
    call assert_small(trim(name), err/maxval(abs(hz)), 1.0e-8_wp)

    call hom%evaluate(z, t + step, hp, dz, dt)
    call hom%evaluate(z, t - step, hm, dz, dt)
    err = norm_inf((hp - hm)/(2.0_wp*step) - ht)
    write(name, "(a,f5.2,sp,f6.2,a)") "dH/dt vs FD, t = ", t, "i"
    call assert_small(trim(name), err/max(1.0_wp, norm_inf(ht)),  &
      1.0e-8_wp)

  end subroutine check_derivatives

end program test_homotopy
