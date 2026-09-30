module polyhomocont__homotopy
!! Homotopies H(z, t) connecting a start system (t = 1) with the target
!! system (t = 0), in projective coordinates on a random affine patch.
!!
!! The path tracker only uses the abstract type `homotopy`: a square
!! system H(z, t) = 0 in N unknowns together with its derivatives with
!! respect to z and t. The concrete `total_degree_homotopy` implements
!!
!!   H_i(z, t) = gamma t G_i(z) + (1 - t) F^h_i(z),   i = 1, ..., n,
!!   H_{n+1}(z) = a . z - 1,
!!
!! with z = (z_1, ..., z_n, w) (homogenizing coordinate last), the
!! homogenized target system F^h, the start system
!! G_i = z_i^{d_i} - b_i w^{d_i} (d_i = degree of F_i) and random complex
!! constants gamma (|gamma| = 1, "gamma trick"), b_i (|b_i| = 1) and a
!! (patch). The random b_i avoid start solutions that coincide with
!! solutions of the target system (e.g. roots of unity), which would
!! give non-generic paths that are constant in t. For generic gamma, the d_1 * ... * d_n
!! paths starting at the solutions of G are smooth for t in (0, 1] and
!! their endpoints at t = 0 include all isolated solutions of F^h = 0
!! (Li, Acta Numerica 6 (1997) 399; HOM4PS-2.0, Lee et al. 2008).
!! Tracking in projective space on the patch a . z = 1 keeps paths that
!! diverge in affine coordinates (solutions at infinity) bounded.
!!
!! The homotopy parameter t is complex: the paths are tracked along the
!! real interval (0, 1], but the end game near t = 0 follows loops
!! |t| = r in the complex t-plane.

  use, intrinsic :: iso_fortran_env, only : int64
  use polyhomocont__config, only : wp
  use polyhomocont__config, only : pi
  use polyhomocont__config, only : status_ok
  use polyhomocont__config, only : err_invalid_input
  use polyhomocont__random, only : rng
  use polyhomocont__system, only : poly_system

  implicit none

  private

  public :: homotopy
  public :: total_degree_homotopy

  type, abstract :: homotopy
    !! A square homotopy H(z, t) = 0, z in C^N, tracked from t = 1 to
    !! t = 0.
  contains
    procedure(hom_nunknowns), deferred :: nunknowns
    procedure(hom_evaluate), deferred :: evaluate
  end type homotopy

  abstract interface

    function hom_nunknowns(self) result(nn)
      !! Number N of unknowns (and equations).
      import :: homotopy
      implicit none
      class(homotopy), intent(in) :: self
      integer :: nn
    end function hom_nunknowns

    subroutine hom_evaluate(self, z, t, h, hz, ht)
      !! Evaluates h = H(z, t), hz = dH/dz (N x N) and ht = dH/dt at a
      !! complex value of t (H is holomorphic in t).
      import :: homotopy, wp
      implicit none
      class(homotopy), intent(inout) :: self
      complex(wp), intent(in) :: z(:)
      complex(wp), intent(in) :: t
      complex(wp), intent(out) :: h(:)
      complex(wp), intent(out) :: hz(:, :)
      complex(wp), intent(out) :: ht(:)
    end subroutine hom_evaluate

  end interface

  type, extends(homotopy) :: total_degree_homotopy
    !! Total-degree (linear) homotopy with the gamma trick on a random
    !! affine patch; see the module documentation.
    class(poly_system), pointer :: sys => null()
      !! Target system (not owned).
    integer :: n = 0
      !! Number of affine variables; the homotopy has n + 1 unknowns.
    integer, allocatable :: deg(:)
      !! Degrees d_i of the target equations.
    complex(wp) :: gamma = (1.0_wp, 0.0_wp)
      !! Random constant of the gamma trick, |gamma| = 1.
    complex(wp), allocatable :: b(:)
      !! Random constants b_i of the start system, |b_i| = 1.
    complex(wp), allocatable :: patch(:)
      !! Coefficients a of the affine patch a . z = 1, unit 2-norm.
    integer :: nstart = 0
      !! Number of start solutions (Bezout number d_1 * ... * d_n).
  contains
    procedure :: init => tdh_init
    procedure :: start_solution => tdh_start_solution
    procedure :: nunknowns => tdh_nunknowns
    procedure :: evaluate => tdh_evaluate
  end type total_degree_homotopy

contains

  subroutine tdh_init(self, sys, gen, status)
    !! Sets up the homotopy for the target system `sys`, drawing gamma
    !! and the patch from `gen`. The homotopy keeps a pointer to `sys`,
    !! which must therefore outlive it. `status` is `err_invalid_input`
    !! if the system is not square (number of degrees /= nvars), has an
    !! equation of degree < 1, or its Bezout number exceeds the default
    !! integer range.

    class(total_degree_homotopy), intent(inout) :: self
    class(poly_system), intent(inout), target :: sys
    type(rng), intent(inout) :: gen
    integer, intent(out) :: status

    integer(int64) :: nstart
    integer :: i

    status = err_invalid_input
    self%sys => sys
    self%n = sys%nvars()
    self%deg = sys%degrees()
    if (self%n < 1) return
    if (size(self%deg) /= self%n) return
    if (any(self%deg < 1)) return

    nstart = 1_int64
    do i = 1, self%n
      nstart = nstart*int(self%deg(i), int64)
      if (nstart > int(huge(1), int64)) return
    end do
    self%nstart = int(nstart)

    self%gamma = gen%unit_complex()
    if (allocated(self%b)) deallocate(self%b)
    allocate(self%b(self%n))
    do i = 1, self%n
      self%b(i) = gen%unit_complex()
    end do
    if (allocated(self%patch)) deallocate(self%patch)
    allocate(self%patch(self%n + 1))
    do i = 1, self%n + 1
      self%patch(i) = gen%normal_complex()
    end do
    self%patch = self%patch/sqrt(sum(abs(self%patch)**2))
    status = status_ok

  end subroutine tdh_init

  subroutine tdh_start_solution(self, k, z)
    !! The k-th start solution, k = 1, ..., nstart: with the mixed-radix
    !! digits j_i of k - 1 (base d_i), w = 1 and
    !! z_i = exp(i (arg b_i + 2 pi j_i) / d_i), a d_i-th root of b_i,
    !! rescaled onto the patch a . z = 1.

    class(total_degree_homotopy), intent(in) :: self
    integer, intent(in) :: k
    complex(wp), intent(out) :: z(:)

    real(wp) :: phi
    integer :: rest
    integer :: digit
    integer :: i

    rest = k - 1
    do i = 1, self%n
      digit = mod(rest, self%deg(i))
      rest = rest/self%deg(i)
      phi = (atan2(aimag(self%b(i)), real(self%b(i), wp))  &
        + 2.0_wp*pi*real(digit, wp))/real(self%deg(i), wp)
      z(i) = cmplx(cos(phi), sin(phi), kind=wp)
    end do
    z(self%n + 1) = (1.0_wp, 0.0_wp)
    z = z/sum(self%patch*z)

  end subroutine tdh_start_solution

  function tdh_nunknowns(self) result(nn)

    class(total_degree_homotopy), intent(in) :: self
    integer :: nn

    nn = self%n + 1

  end function tdh_nunknowns

  subroutine tdh_evaluate(self, z, t, h, hz, ht)

    class(total_degree_homotopy), intent(inout) :: self
    complex(wp), intent(in) :: z(:)
    complex(wp), intent(in) :: t
    complex(wp), intent(out) :: h(:)
    complex(wp), intent(out) :: hz(:, :)
    complex(wp), intent(out) :: ht(:)

    complex(wp) :: f(self%n)
    complex(wp) :: jf(self%n, self%n + 1)
    complex(wp) :: g
    complex(wp) :: gt
    complex(wp) :: omt
    complex(wp) :: zd1
    complex(wp) :: wd1
    complex(wp) :: w
    integer :: n
    integer :: d
    integer :: i

    n = self%n
    w = z(n + 1)
    call self%sys%evaluate_homogeneous(z, f, jf)

    gt = self%gamma*t
    omt = 1.0_wp - t
    hz(1:n, :) = omt*jf
    do i = 1, n
      d = self%deg(i)
      zd1 = z(i)**(d - 1)
      wd1 = w**(d - 1)
      g = zd1*z(i) - self%b(i)*wd1*w
      h(i) = gt*g + omt*f(i)
      ht(i) = self%gamma*g - f(i)
      hz(i, i) = hz(i, i) + gt*d*zd1
      hz(i, n + 1) = hz(i, n + 1) - gt*self%b(i)*d*wd1
    end do

    h(n + 1) = sum(self%patch*z) - 1.0_wp
    hz(n + 1, :) = self%patch
    ht(n + 1) = (0.0_wp, 0.0_wp)

  end subroutine tdh_evaluate

end module polyhomocont__homotopy
