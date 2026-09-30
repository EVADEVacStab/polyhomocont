module polyhomocont__linalg
!! Dense complex linear algebra for the (small) Newton and predictor
!! systems of the path tracker: LU factorization with partial pivoting,
!! forward/back substitution and the 1-norm condition number.
!!
!! The systems solved by polyhomocont have dimension n + 1, with n the
!! number of variables, which is small in practice; a self-contained
!! implementation is therefore used instead of LAPACK.

  use polyhomocont__config, only : wp
  use polyhomocont__config, only : status_ok
  use polyhomocont__config, only : err_singular_matrix

  implicit none

  private

  public :: lu_factor
  public :: lu_solve
  public :: solve_linear
  public :: cond1
  public :: norm_inf

contains

  subroutine lu_factor(a, ipiv, status)
    !! In-place LU factorization P*A = L*U with partial pivoting (the
    !! pivot is the entry of largest modulus in the column). On exit, the
    !! strict lower triangle of `a` holds L (unit diagonal implied) and the
    !! upper triangle holds U; row i was interchanged with row ipiv(i) in
    !! step i. `status` is `err_singular_matrix` if a pivot is exactly
    !! zero.

    complex(wp), intent(inout) :: a(:, :)
    integer, intent(out) :: ipiv(:)
    integer, intent(out) :: status

    complex(wp) :: tmp(size(a, 2))
    complex(wp) :: piv
    real(wp) :: amax
    integer :: n
    integer :: i
    integer :: k
    integer :: p

    status = status_ok
    n = size(a, 1)

    do k = 1, n
      ! Find the pivot row
      p = k
      amax = abs(a(k, k))
      do i = k + 1, n
        if (abs(a(i, k)) > amax) then
          amax = abs(a(i, k))
          p = i
        end if
      end do
      ipiv(k) = p
      if (amax == 0.0_wp) then
        status = err_singular_matrix
        return
      end if
      if (p /= k) then
        tmp = a(k, :)
        a(k, :) = a(p, :)
        a(p, :) = tmp
      end if

      ! Eliminate below the pivot
      piv = a(k, k)
      do i = k + 1, n
        a(i, k) = a(i, k)/piv
        a(i, k+1:n) = a(i, k+1:n) - a(i, k)*a(k, k+1:n)
      end do
    end do

  end subroutine lu_factor

  subroutine lu_solve(a, ipiv, b)
    !! Solves A*x = b given the factorization from `lu_factor`. On exit,
    !! `b` holds the solution x.

    complex(wp), intent(in) :: a(:, :)
    integer, intent(in) :: ipiv(:)
    complex(wp), intent(inout) :: b(:)

    complex(wp) :: tmp
    integer :: n
    integer :: i
    integer :: k

    n = size(a, 1)

    ! Apply the row interchanges
    do k = 1, n
      if (ipiv(k) /= k) then
        tmp = b(k)
        b(k) = b(ipiv(k))
        b(ipiv(k)) = tmp
      end if
    end do

    ! Forward substitution with L (unit diagonal)
    do i = 2, n
      b(i) = b(i) - sum(a(i, 1:i-1)*b(1:i-1))
    end do

    ! Back substitution with U
    do i = n, 1, -1
      b(i) = (b(i) - sum(a(i, i+1:n)*b(i+1:n)))/a(i, i)
    end do

  end subroutine lu_solve

  subroutine solve_linear(a, b, x, status)
    !! Solves A*x = b without modifying `a` or `b`.

    complex(wp), intent(in) :: a(:, :)
    complex(wp), intent(in) :: b(:)
    complex(wp), intent(out) :: x(:)
    integer, intent(out) :: status

    complex(wp) :: lu(size(a, 1), size(a, 2))
    integer :: ipiv(size(a, 1))

    lu = a
    call lu_factor(lu, ipiv, status)
    if (status /= status_ok) return
    x = b
    call lu_solve(lu, ipiv, x)

  end subroutine solve_linear

  function cond1(a) result(c)
    !! Condition number of `a` in the 1-norm, ||A||_1 ||A^-1||_1, with the
    !! inverse computed explicitly (the matrices are small). Returns
    !! `huge(1.0_wp)` for an exactly singular matrix.

    complex(wp), intent(in) :: a(:, :)
    real(wp) :: c

    complex(wp) :: lu(size(a, 1), size(a, 2))
    complex(wp) :: col(size(a, 1))
    integer :: ipiv(size(a, 1))
    real(wp) :: anorm
    real(wp) :: ainvnorm
    integer :: status
    integer :: n
    integer :: j

    n = size(a, 1)
    anorm = 0.0_wp
    do j = 1, n
      anorm = max(anorm, sum(abs(a(:, j))))
    end do

    lu = a
    call lu_factor(lu, ipiv, status)
    if (status /= status_ok) then
      c = huge(1.0_wp)
      return
    end if

    ainvnorm = 0.0_wp
    do j = 1, n
      col = (0.0_wp, 0.0_wp)
      col(j) = (1.0_wp, 0.0_wp)
      call lu_solve(lu, ipiv, col)
      ainvnorm = max(ainvnorm, sum(abs(col)))
    end do
    c = anorm*ainvnorm

  end function cond1

  pure function norm_inf(z) result(r)
    !! Maximum modulus of the components of `z` (0 for an empty vector).

    complex(wp), intent(in) :: z(:)
    real(wp) :: r

    integer :: i

    r = 0.0_wp
    do i = 1, size(z)
      r = max(r, abs(z(i)))
    end do

  end function norm_inf

end module polyhomocont__linalg
