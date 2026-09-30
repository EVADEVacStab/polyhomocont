module polyhomocont__sparse
!! Built-in sparse representation of polynomial systems.
!!
!! Each equation is a list of terms c * x^a with complex coefficients c
!! and exponent vectors a (columns of an integer matrix). Equations are
!! added one at a time:
!!
!!   type(sparse_system) :: sys
!!   call sys%init(2)
!!   ! x1^2 + x2^2 - 5
!!   call sys%add_equation([1.0_wp, 1.0_wp, -5.0_wp],  &
!!     reshape([2, 0,  0, 2,  0, 0], [2, 3]), status)
!!
!! F and dF/dx are evaluated together from a table of powers of the
!! variables; the derivative of each term is formed from prefix and
!! suffix products, so no division by x_j is needed and points with
!! zero coordinates are handled exactly. The homogenized system is
!! evaluated exactly as a polynomial in (z, w), also at infinity (w = 0).

  use polyhomocont__config, only : wp
  use polyhomocont__config, only : status_ok
  use polyhomocont__config, only : err_invalid_input
  use polyhomocont__system, only : poly_system

  implicit none

  private

  public :: sparse_system

  type, extends(poly_system) :: sparse_system
    !! Sparse polynomial system. The terms of equation i are
    !! first(i):first(i+1)-1 in `coefs`, `exps` and `hexps`.
    private
    integer :: n = 0
      !! Number of variables.
    integer :: neq = 0
      !! Number of equations added so far.
    integer :: maxexp = 0
      !! Largest exponent of any variable (including w in `hexps`).
    integer, allocatable :: first(:)
    integer, allocatable :: deg(:)
    integer, allocatable :: exps(:, :)
      !! Exponent vectors, shape (n, nterms).
    integer, allocatable :: hexps(:, :)
      !! Exponent vectors of the homogenized terms, shape (n + 1, nterms);
      !! the last row is the exponent d_i - |a| of w.
    complex(wp), allocatable :: coefs(:)
  contains
    procedure :: init => sparse_init
    procedure, private :: sparse_add_equation_complex
    procedure, private :: sparse_add_equation_real
    generic :: add_equation => sparse_add_equation_complex,  &
      sparse_add_equation_real
    procedure :: nequations => sparse_nequations
    procedure :: nterms => sparse_nterms
    procedure :: nvars => sparse_nvars
    procedure :: degrees => sparse_degrees
    procedure :: evaluate => sparse_evaluate
    procedure :: evaluate_homogeneous => sparse_evaluate_homogeneous
  end type sparse_system

contains

  subroutine sparse_init(self, n)
    !! Initializes an empty system in n variables; any previous content is
    !! discarded. Afterwards, add n equations with `add_equation`.

    class(sparse_system), intent(inout) :: self
    integer, intent(in) :: n

    self%n = n
    self%neq = 0
    self%maxexp = 0
    self%first = [1]
    self%deg = [integer ::]
    self%exps = reshape([integer ::], [n, 0])
    self%hexps = reshape([integer ::], [n + 1, 0])
    self%coefs = [complex(wp) ::]

  end subroutine sparse_init

  subroutine sparse_add_equation_complex(self, coefs, exps, status)
    !! Appends the equation sum_k coefs(k) * x^exps(:, k) = 0. Terms with
    !! zero coefficient are dropped; repeated exponent vectors are allowed
    !! (their coefficients effectively add up). `status` is
    !! `err_invalid_input` if the shapes do not match, an exponent is
    !! negative, the system already has n equations, or the equation is
    !! constant (degree 0).

    class(sparse_system), intent(inout) :: self
    complex(wp), intent(in) :: coefs(:)
    integer, intent(in) :: exps(:, :)
    integer, intent(out) :: status

    logical :: keep(size(coefs))
    integer, allocatable :: newexps(:, :)
    integer, allocatable :: newhexps(:, :)
    integer :: nkeep
    integer :: d
    integer :: k
    integer :: m

    status = err_invalid_input
    if (.not. allocated(self%first)) return
    if (self%neq >= self%n) return
    if (size(exps, 1) /= self%n) return
    if (size(exps, 2) /= size(coefs)) return
    if (any(exps < 0)) return

    keep = coefs /= (0.0_wp, 0.0_wp)
    nkeep = count(keep)
    d = 0
    do k = 1, size(coefs)
      if (keep(k)) d = max(d, sum(exps(:, k)))
    end do
    if (d < 1) return

    ! Append the kept terms
    allocate(newexps(self%n, nkeep))
    allocate(newhexps(self%n + 1, nkeep))
    m = 0
    do k = 1, size(coefs)
      if (.not. keep(k)) cycle
      m = m + 1
      newexps(:, m) = exps(:, k)
      newhexps(1:self%n, m) = exps(:, k)
      newhexps(self%n + 1, m) = d - sum(exps(:, k))
    end do
    self%exps = reshape([self%exps, newexps],  &
      [self%n, size(self%exps, 2) + nkeep])
    self%hexps = reshape([self%hexps, newhexps],  &
      [self%n + 1, size(self%hexps, 2) + nkeep])
    self%coefs = [self%coefs, pack(coefs, keep)]
    self%neq = self%neq + 1
    self%first = [self%first, self%first(self%neq) + nkeep]
    self%deg = [self%deg, d]
    self%maxexp = max(self%maxexp, maxval(newhexps))
    status = status_ok

  end subroutine sparse_add_equation_complex

  subroutine sparse_add_equation_real(self, coefs, exps, status)
    !! Same as the complex version, for real coefficients.

    class(sparse_system), intent(inout) :: self
    real(wp), intent(in) :: coefs(:)
    integer, intent(in) :: exps(:, :)
    integer, intent(out) :: status

    call self%sparse_add_equation_complex(  &
      cmplx(coefs, 0.0_wp, kind=wp), exps, status)

  end subroutine sparse_add_equation_real

  function sparse_nequations(self) result(neq)
    !! Number of equations added so far.

    class(sparse_system), intent(in) :: self
    integer :: neq

    neq = self%neq

  end function sparse_nequations

  function sparse_nterms(self, i) result(nt)
    !! Number of (nonzero) terms of equation i.

    class(sparse_system), intent(in) :: self
    integer, intent(in) :: i
    integer :: nt

    nt = self%first(i+1) - self%first(i)

  end function sparse_nterms

  function sparse_nvars(self) result(n)

    class(sparse_system), intent(in) :: self
    integer :: n

    n = self%n

  end function sparse_nvars

  function sparse_degrees(self) result(d)
    !! Degrees of the equations added so far (the solver checks that
    !! there are n of them).

    class(sparse_system), intent(in) :: self
    integer, allocatable :: d(:)

    d = self%deg

  end function sparse_degrees

  subroutine sparse_evaluate(self, x, f, jac)

    class(sparse_system), intent(inout) :: self
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(out) :: f(:)
    complex(wp), intent(out) :: jac(:, :)

    call eval_terms(self%neq, self%first, self%exps, self%coefs,  &
      self%maxexp, x, f, jac)

  end subroutine sparse_evaluate

  subroutine sparse_evaluate_homogeneous(self, z, f, jac)
    !! Exact evaluation of the homogenized system and its Jacobian at
    !! z = (z_1, ..., z_n, w), valid also for w = 0.

    class(sparse_system), intent(inout) :: self
    complex(wp), intent(in) :: z(:)
    complex(wp), intent(out) :: f(:)
    complex(wp), intent(out) :: jac(:, :)

    call eval_terms(self%neq, self%first, self%hexps, self%coefs,  &
      self%maxexp, z, f, jac)

  end subroutine sparse_evaluate_homogeneous

  subroutine eval_terms(neq, first, exps, coefs, maxexp, v, f, jac)
    !! Evaluates f(i) = sum_k c_k v^e_k over the terms of equation i and
    !! jac(i, j) = df_i/dv_j, for m = size(v) variables. For each term,
    !! with p_j = v_j^e_j, the derivative with respect to v_j is
    !! c e_j v_j^(e_j - 1) prod_{l /= j} p_l, where the product over
    !! l /= j is (prefix product up to j - 1) * (suffix product from j + 1).

    integer, intent(in) :: neq
    integer, intent(in) :: first(:)
    integer, intent(in) :: exps(:, :)
    complex(wp), intent(in) :: coefs(:)
    integer, intent(in) :: maxexp
    complex(wp), intent(in) :: v(:)
    complex(wp), intent(out) :: f(:)
    complex(wp), intent(out) :: jac(:, :)

    complex(wp) :: pw(0:maxexp, size(v))
    complex(wp) :: pre(0:size(v))
    complex(wp) :: suf(1:size(v)+1)
    complex(wp) :: c
    integer :: m
    integer :: i
    integer :: j
    integer :: k
    integer :: e

    m = size(v)

    ! Power table pw(p, j) = v_j^p
    pw(0, :) = (1.0_wp, 0.0_wp)
    do e = 1, maxexp
      pw(e, :) = pw(e-1, :)*v
    end do

    f = (0.0_wp, 0.0_wp)
    jac = (0.0_wp, 0.0_wp)
    do i = 1, neq
      do k = first(i), first(i+1) - 1
        c = coefs(k)
        pre(0) = (1.0_wp, 0.0_wp)
        do j = 1, m
          pre(j) = pre(j-1)*pw(exps(j, k), j)
        end do
        suf(m+1) = (1.0_wp, 0.0_wp)
        do j = m, 1, -1
          suf(j) = pw(exps(j, k), j)*suf(j+1)
        end do
        f(i) = f(i) + c*pre(m)
        do j = 1, m
          e = exps(j, k)
          if (e > 0) then
            jac(i, j) = jac(i, j) + c*e*pw(e-1, j)*pre(j-1)*suf(j+1)
          end if
        end do
      end do
    end do

  end subroutine eval_terms

end module polyhomocont__sparse
