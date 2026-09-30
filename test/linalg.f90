program test_linalg
!! Tests of the complex LU factorization, linear solves and condition
!! numbers.

  use polyhomocont, only : wp
  use polyhomocont, only : err_singular_matrix
  use polyhomocont__random, only : rng
  use polyhomocont__linalg, only : lu_factor
  use polyhomocont__linalg, only : lu_solve
  use polyhomocont__linalg, only : solve_linear
  use polyhomocont__linalg, only : cond1
  use polyhomocont__linalg, only : norm_inf
  use testing, only : assert_status
  use testing, only : assert_small
  use testing, only : assert_close
  use testing, only : assert_true
  use testing, only : finish

  implicit none

  type(rng) :: gen
  integer :: n
  character(len=40) :: name

  call gen%seed(2026)

  ! Random dense systems of various sizes
  do n = 1, 20, 3
    write(name, "(a,i0)") "random solve, n = ", n
    call check_random(n, trim(name))
  end do

  call check_pivoting()
  call check_singular()
  call check_cond()

  call finish()

contains

  subroutine check_random(n, name)
    !! Builds b = A*xref for random A and xref, solves A*x = b and checks
    !! the error relative to xref.

    integer, intent(in) :: n
    character(len=*), intent(in) :: name

    complex(wp) :: a(n, n)
    complex(wp) :: xref(n)
    complex(wp) :: b(n)
    complex(wp) :: x(n)
    integer :: status
    integer :: i
    integer :: j

    do j = 1, n
      do i = 1, n
        a(i, j) = gen%normal_complex()
      end do
      xref(j) = gen%normal_complex()
    end do
    b = matmul(a, xref)

    call solve_linear(a, b, x, status)
    call assert_status(name//": status", status)
    call assert_small(name//": error",  &
      norm_inf(x - xref)/norm_inf(xref), 1.0e-11_wp)

  end subroutine check_random

  subroutine check_pivoting()
    !! A matrix with a zero in the (1, 1) position requires a row
    !! interchange; reusing the factorization for two right-hand sides.

    complex(wp) :: a(3, 3)
    complex(wp) :: lu(3, 3)
    complex(wp) :: xref(3)
    complex(wp) :: b(3)
    integer :: ipiv(3)
    integer :: status

    a = reshape([  &
      (0.0_wp, 0.0_wp), (2.0_wp, 1.0_wp), (1.0_wp, 0.0_wp),  &
      (1.0_wp, -1.0_wp), (0.0_wp, 0.0_wp), (3.0_wp, 0.0_wp),  &
      (4.0_wp, 0.0_wp), (1.0_wp, 1.0_wp), (0.0_wp, 2.0_wp)], [3, 3])
    lu = a
    call lu_factor(lu, ipiv, status)
    call assert_status("pivoting: factor", status)

    xref = [(1.0_wp, 0.0_wp), (0.0_wp, 1.0_wp), (-2.0_wp, 0.5_wp)]
    b = matmul(a, xref)
    call lu_solve(lu, ipiv, b)
    call assert_small("pivoting: error rhs 1", norm_inf(b - xref),  &
      1.0e-14_wp)

    xref = [(0.5_wp, 0.5_wp), (3.0_wp, 0.0_wp), (0.0_wp, -1.0_wp)]
    b = matmul(a, xref)
    call lu_solve(lu, ipiv, b)
    call assert_small("pivoting: error rhs 2", norm_inf(b - xref),  &
      1.0e-14_wp)

  end subroutine check_pivoting

  subroutine check_singular()
    !! An exactly singular matrix (two equal rows) must be reported, and
    !! its condition number must be huge.

    complex(wp) :: a(3, 3)
    complex(wp) :: b(3)
    complex(wp) :: x(3)
    integer :: status

    a = reshape([  &
      (1.0_wp, 0.0_wp), (1.0_wp, 0.0_wp), (0.0_wp, 1.0_wp),  &
      (2.0_wp, 0.0_wp), (2.0_wp, 0.0_wp), (1.0_wp, 0.0_wp),  &
      (0.0_wp, 3.0_wp), (0.0_wp, 3.0_wp), (5.0_wp, 0.0_wp)], [3, 3])
    b = (1.0_wp, 0.0_wp)
    call solve_linear(a, b, x, status)
    call assert_true("singular: detected", status == err_singular_matrix)
    call assert_true("singular: cond1 = huge", cond1(a) == huge(1.0_wp))

  end subroutine check_singular

  subroutine check_cond()
    !! Condition numbers with known values.

    complex(wp) :: d(3, 3)
    complex(wp) :: a(2, 2)
    real(wp) :: eps

    ! Diagonal matrix: cond1 = max|d_i| / min|d_i|
    d = (0.0_wp, 0.0_wp)
    d(1, 1) = (1.0e3_wp, 0.0_wp)
    d(2, 2) = (0.0_wp, -2.0_wp)
    d(3, 3) = (0.3_wp, 0.4_wp)
    call assert_close("cond1: diagonal", cond1(d), 1.0e3_wp/0.5_wp,  &
      1.0e-13_wp)

    ! A = [[1, 1], [1, 1 + eps]]: ||A||_1 = 2 + eps,
    ! A^-1 = [[1 + eps, -1], [-1, 1]]/eps, ||A^-1||_1 = (2 + eps)/eps
    eps = 1.0e-6_wp
    a = reshape([(1.0_wp, 0.0_wp), (1.0_wp, 0.0_wp),  &
      (1.0_wp, 0.0_wp), cmplx(1.0_wp + eps, 0.0_wp, kind=wp)], [2, 2])
    call assert_close("cond1: nearly singular", cond1(a),  &
      (2.0_wp + eps)**2/eps, 1.0e-8_wp)

  end subroutine check_cond

end program test_linalg
