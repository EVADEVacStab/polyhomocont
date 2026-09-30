module polyhomocont__system
!! Abstract interface for square polynomial systems F(x) = 0,
!! F: C^n -> C^n.
!!
!! Users extend `poly_system` and implement
!!   - `nvars`:    the number of variables n (= number of equations),
!!   - `degrees`:  the total degree d_i of each equation,
!!   - `evaluate`: F(x) together with the Jacobian dF/dx.
!! The built-in `sparse_system` (module polyhomocont__sparse) is one such
!! extension; systems given in closed form can implement their own.
!!
!! The path tracker works in projective coordinates z = (z_1, ..., z_n, w)
!! with the homogenizing coordinate w last, x = z(1:n)/w, and evaluates
!! the homogenized system F^h_i(z) = w^{d_i} F_i(z(1:n)/w) through
!! `evaluate_homogeneous`. Its default implementation below maps to the
!! affine evaluator; it is exact for w /= 0 but loses accuracy as
!! w -> 0 (solutions at infinity), since x then becomes large.
!! Extensions that can evaluate F^h directly (like `sparse_system`)
!! should override it.

  use polyhomocont__config, only : wp
  use polyhomocont__linalg, only : norm_inf

  implicit none

  private

  public :: poly_system

  type, abstract :: poly_system
    !! A square system of n polynomial equations in n complex variables.
  contains
    procedure(sys_nvars), deferred :: nvars
    procedure(sys_degrees), deferred :: degrees
    procedure(sys_evaluate), deferred :: evaluate
    procedure :: evaluate_homogeneous => poly_system_evaluate_homogeneous
  end type poly_system

  abstract interface

    function sys_nvars(self) result(n)
      !! Number of variables (and equations) n.
      import :: poly_system
      implicit none
      class(poly_system), intent(in) :: self
      integer :: n
    end function sys_nvars

    function sys_degrees(self) result(d)
      !! Total degrees d_i >= 1 of the n equations. They must be exact:
      !! the start system of the total-degree homotopy and the
      !! homogenization are built from them.
      import :: poly_system
      implicit none
      class(poly_system), intent(in) :: self
      integer, allocatable :: d(:)
    end function sys_degrees

    subroutine sys_evaluate(self, x, f, jac)
      !! Evaluates f = F(x) and jac(i, j) = dF_i/dx_j at x(1:n).
      import :: poly_system, wp
      implicit none
      class(poly_system), intent(inout) :: self
      complex(wp), intent(in) :: x(:)
      complex(wp), intent(out) :: f(:)
      complex(wp), intent(out) :: jac(:, :)
    end subroutine sys_evaluate

  end interface

contains

  subroutine poly_system_evaluate_homogeneous(self, z, f, jac)
    !! Evaluates the homogenized system f = F^h(z) and its Jacobian
    !! jac(i, j) = dF^h_i/dz_j (shape n x (n + 1)) at z(1:n+1), with the
    !! homogenizing coordinate w = z(n+1) last. With x = z(1:n)/w:
    !!   F^h_i        = w^{d_i} F_i(x)
    !!   dF^h_i/dz_j  = w^{d_i - 1} dF_i/dx_j                (j <= n)
    !!   dF^h_i/dw    = w^{d_i - 1} (d_i F_i(x) - sum_j x_j dF_i/dx_j)
    !! Exactly w = 0 is replaced by a tiny w (relative size epsilon), which
    !! only happens for non-generic input.

    class(poly_system), intent(inout) :: self
    complex(wp), intent(in) :: z(:)
    complex(wp), intent(out) :: f(:)
    complex(wp), intent(out) :: jac(:, :)

    integer, allocatable :: d(:)
    complex(wp), allocatable :: x(:)
    complex(wp), allocatable :: fa(:)
    complex(wp), allocatable :: ja(:, :)
    complex(wp) :: w
    complex(wp) :: wd1
    integer :: n
    integer :: i

    n = size(z) - 1
    d = self%degrees()
    allocate(fa(n))
    allocate(ja(n, n))

    w = z(n+1)
    if (w == (0.0_wp, 0.0_wp)) then
      w = cmplx(epsilon(1.0_wp)*max(norm_inf(z), tiny(1.0_wp)), 0.0_wp,  &
        kind=wp)
    end if
    x = z(1:n)/w
    call self%evaluate(x, fa, ja)

    do i = 1, n
      wd1 = w**(d(i) - 1)
      f(i) = wd1*w*fa(i)
      jac(i, 1:n) = wd1*ja(i, :)
      jac(i, n+1) = wd1*(d(i)*fa(i) - sum(x*ja(i, :)))
    end do

  end subroutine poly_system_evaluate_homogeneous

end module polyhomocont__system
