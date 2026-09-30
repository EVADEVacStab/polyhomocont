module testing
!! Assertion helpers shared by the test programs. Each helper prints
!! `ok   <name>` or `FAIL <name>` and counts the failures; a test program
!! finishes with `call finish()`, which stops with a nonzero exit code if
!! any assertion failed.

  use polyhomocont, only : wp
  use polyhomocont, only : status_ok
  use polyhomocont, only : poly_system

  implicit none

  private

  public :: assert_status
  public :: assert_close
  public :: assert_small
  public :: assert_true
  public :: assert_equal
  public :: finish
  public :: jacobian_fd_error
  public :: homogeneous_jacobian_fd_error
  public :: homogeneity_error

  integer :: nfail = 0
  integer :: npass = 0

contains

  subroutine assert_status(name, status)

    character(len=*), intent(in) :: name
    integer, intent(in) :: status

    if (status /= status_ok) then
      print "(a,a,a,i0)", "FAIL ", name, ": status = ", status
      nfail = nfail + 1
    else
      npass = npass + 1
    end if

  end subroutine assert_status

  subroutine assert_close(name, val, ref, rtol)
    !! Checks |val - ref| <= rtol*|ref|.

    character(len=*), intent(in) :: name
    real(wp), intent(in) :: val
    real(wp), intent(in) :: ref
    real(wp), intent(in) :: rtol

    if (abs(val - ref) > rtol*abs(ref)) then
      print "(a,a,a,es16.8,a,es16.8,a,es9.2)", "FAIL ", name, ": got ",  &
        val, ", expected ", ref, ", rtol ", rtol
      nfail = nfail + 1
    else
      print "(a,a,a,es16.8,a,es16.8,a)", "ok   ", name, ": got ", val,  &
        " (ref ", ref, ")"
      npass = npass + 1
    end if

  end subroutine assert_close

  subroutine assert_small(name, val, tol)
    !! Checks |val| <= tol, e.g. for residuals and errors.

    character(len=*), intent(in) :: name
    real(wp), intent(in) :: val
    real(wp), intent(in) :: tol

    if (.not. (abs(val) <= tol)) then
      print "(a,a,a,es10.3,a,es9.2)", "FAIL ", name, ": ", val,  &
        " > tol ", tol
      nfail = nfail + 1
    else
      print "(a,a,a,es10.3,a,es9.2)", "ok   ", name, ": ", val,  &
        " <= ", tol
      npass = npass + 1
    end if

  end subroutine assert_small

  subroutine assert_true(name, cond)

    character(len=*), intent(in) :: name
    logical, intent(in) :: cond

    if (.not. cond) then
      print "(a,a)", "FAIL ", name
      nfail = nfail + 1
    else
      print "(a,a)", "ok   ", name
      npass = npass + 1
    end if

  end subroutine assert_true

  subroutine assert_equal(name, val, ref)

    character(len=*), intent(in) :: name
    integer, intent(in) :: val
    integer, intent(in) :: ref

    if (val /= ref) then
      print "(a,a,a,i0,a,i0)", "FAIL ", name, ": got ", val,  &
        ", expected ", ref
      nfail = nfail + 1
    else
      print "(a,a,a,i0)", "ok   ", name, ": ", val
      npass = npass + 1
    end if

  end subroutine assert_equal

  subroutine finish()
    !! Prints a summary and stops with exit code 1 if anything failed.

    print "(i0,a,i0,a)", npass, " passed, ", nfail, " failed"
    if (nfail > 0) error stop 1

  end subroutine finish

  function jacobian_fd_error(sys, x) result(err)
    !! Maximum deviation between the Jacobian returned by sys%evaluate and
    !! central finite differences (polynomials are holomorphic, so a real
    !! step in each complex variable suffices), relative to the largest
    !! Jacobian entry.

    class(poly_system), intent(inout) :: sys
    complex(wp), intent(in) :: x(:)
    real(wp) :: err

    complex(wp) :: f(size(x))
    complex(wp) :: fp(size(x))
    complex(wp) :: fm(size(x))
    complex(wp) :: jac(size(x), size(x))
    complex(wp) :: jdum(size(x), size(x))
    complex(wp) :: xp(size(x))
    real(wp) :: h
    integer :: j

    call sys%evaluate(x, f, jac)
    err = 0.0_wp
    do j = 1, size(x)
      h = 1.0e-5_wp*max(1.0_wp, abs(x(j)))
      xp = x
      xp(j) = x(j) + h
      call sys%evaluate(xp, fp, jdum)
      xp(j) = x(j) - h
      call sys%evaluate(xp, fm, jdum)
      err = max(err, maxval(abs((fp - fm)/(2.0_wp*h) - jac(:, j))))
    end do
    err = err/max(1.0_wp, maxval(abs(jac)))

  end function jacobian_fd_error

  function homogeneous_jacobian_fd_error(sys, z) result(err)
    !! Like jacobian_fd_error, for sys%evaluate_homogeneous at the
    !! projective point z(1:n+1).

    class(poly_system), intent(inout) :: sys
    complex(wp), intent(in) :: z(:)
    real(wp) :: err

    complex(wp) :: f(size(z) - 1)
    complex(wp) :: fp(size(z) - 1)
    complex(wp) :: fm(size(z) - 1)
    complex(wp) :: jac(size(z) - 1, size(z))
    complex(wp) :: jdum(size(z) - 1, size(z))
    complex(wp) :: zp(size(z))
    real(wp) :: h
    integer :: j

    call sys%evaluate_homogeneous(z, f, jac)
    err = 0.0_wp
    do j = 1, size(z)
      h = 1.0e-5_wp*max(1.0_wp, abs(z(j)))
      zp = z
      zp(j) = z(j) + h
      call sys%evaluate_homogeneous(zp, fp, jdum)
      zp(j) = z(j) - h
      call sys%evaluate_homogeneous(zp, fm, jdum)
      err = max(err, maxval(abs((fp - fm)/(2.0_wp*h) - jac(:, j))))
    end do
    err = err/max(1.0_wp, maxval(abs(jac)))

  end function homogeneous_jacobian_fd_error

  function homogeneity_error(sys, z, lambda) result(err)
    !! Checks F^h_i(lambda*z) = lambda^{d_i} F^h_i(z) and
    !! F^h_i(z) = w^{d_i} F_i(z(1:n)/w) (w = z(n+1) /= 0); returns the
    !! largest relative deviation.

    class(poly_system), intent(inout) :: sys
    complex(wp), intent(in) :: z(:)
    complex(wp), intent(in) :: lambda
    real(wp) :: err

    complex(wp) :: f(size(z) - 1)
    complex(wp) :: fl(size(z) - 1)
    complex(wp) :: fa(size(z) - 1)
    complex(wp) :: jh(size(z) - 1, size(z))
    complex(wp) :: ja(size(z) - 1, size(z) - 1)
    integer, allocatable :: d(:)
    integer :: n

    n = size(z) - 1
    d = sys%degrees()
    call sys%evaluate_homogeneous(z, f, jh)
    call sys%evaluate_homogeneous(lambda*z, fl, jh)
    call sys%evaluate(z(1:n)/z(n+1), fa, ja)
    err = maxval(abs(fl - lambda**d*f))/max(1.0_wp, maxval(abs(fl)))
    err = max(err,  &
      maxval(abs(f - z(n+1)**d*fa))/max(1.0_wp, maxval(abs(f))))

  end function homogeneity_error

end module testing
