module testing
!! Assertion helpers shared by the test programs. Each helper prints
!! `ok   <name>` or `FAIL <name>` and counts the failures; a test program
!! finishes with `call finish()`, which stops with a nonzero exit code if
!! any assertion failed.

  use polyhomocont, only : wp
  use polyhomocont, only : status_ok

  implicit none

  private

  public :: assert_status
  public :: assert_close
  public :: assert_small
  public :: assert_true
  public :: assert_equal
  public :: finish

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

end module testing
