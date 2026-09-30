module polyhomocont__solve
!! Driver that computes all isolated solutions of a square polynomial
!! system with the total-degree homotopy, and post-processes the path
!! endpoints.
!!
!! For each of the d_1 * ... * d_n start solutions, the path is tracked
!! from t = 1 to t = 0 in projective coordinates. Each endpoint z is
!! refined with Newton's method on the target system and classified:
!!   - at infinity if |w| <= infinity_tol * ||z|| (w = z(n+1)), otherwise
!!     finite with x = z(1:n)/w;
!!   - singular if the condition number of the Jacobian of the
!!     (patched, homogenized) target system exceeds singular_cond;
!!   - real if ||Im x|| <= real_tol * max(1, ||x||).
!! The distinct finite nonsingular endpoints are collected in
!! `solutions`. Two nonsingular paths cannot end at the same point for a
!! generic homotopy, so such duplicates indicate path jumping and are
!! counted in `nduplicates`.
!!
!! Singular endpoints are only approached with limited accuracy (and the
!! tracker may fail very close to t = 0) until an end game is available.

  use polyhomocont__config, only : wp
  use polyhomocont__config, only : status_ok
  use polyhomocont__random, only : rng
  use polyhomocont__linalg, only : lu_factor
  use polyhomocont__linalg, only : lu_solve
  use polyhomocont__linalg, only : cond1
  use polyhomocont__linalg, only : norm_inf
  use polyhomocont__system, only : poly_system
  use polyhomocont__homotopy, only : total_degree_homotopy
  use polyhomocont__tracker, only : tracker_options
  use polyhomocont__tracker, only : path_info
  use polyhomocont__tracker, only : track_path
  use polyhomocont__tracker, only : path_success

  implicit none

  private

  public :: solve_options
  public :: path_result
  public :: solve_result
  public :: solve

  type :: solve_options
    !! Parameters of the solver.
    type(tracker_options) :: tracker
      !! Parameters of the path tracker.
    integer :: seed = 1
      !! Seed for the random constants (gamma, patch). Results are
      !! reproducible for a fixed seed.
    integer :: refine_maxit = 5
      !! Maximum number of Newton iterations refining each endpoint.
    real(wp) :: infinity_tol = 1.0e-8_wp
      !! Endpoints with |w| <= infinity_tol * ||z|| are at infinity.
    real(wp) :: singular_cond = 1.0e10_wp
      !! Endpoints whose Jacobian has a larger condition number are
      !! singular.
    real(wp) :: real_tol = 1.0e-8_wp
      !! Solutions with ||Im x|| <= real_tol * max(1, ||x||) are real.
    real(wp) :: dedup_tol = 1.0e-8_wp
      !! Relative distance below which two solutions are identical.
  end type solve_options

  type :: path_result
    !! Endpoint and classification of one path.
    integer :: status = path_success
      !! Tracker status (path_success or path_failed_*).
    real(wp) :: t = 1.0_wp
      !! Value of t reached (0 on success).
    complex(wp), allocatable :: z(:)
      !! Endpoint in projective coordinates (z_1, ..., z_n, w).
    complex(wp), allocatable :: x(:)
      !! Affine endpoint z(1:n)/w; allocated only if the endpoint is
      !! finite.
    logical :: at_infinity = .false.
    logical :: singular = .false.
    logical :: is_real = .false.
    real(wp) :: cond = 0.0_wp
      !! 1-norm condition number of the Jacobian at the endpoint.
    real(wp) :: residual = 0.0_wp
      !! max_i |F^h_i(z)| / ||z||^{d_i}.
    integer :: solution_index = 0
      !! Index into solve_result%solutions (0 if the endpoint is not a
      !! finite nonsingular solution).
    integer :: naccepted = 0
    integer :: nrejected = 0
  end type path_result

  type :: solve_result
    !! Result of `solve`.
    integer :: n = 0
      !! Number of variables.
    integer :: npaths = 0
      !! Number of tracked paths (Bezout number).
    type(path_result), allocatable :: paths(:)
    integer :: nsolutions = 0
      !! Number of distinct finite nonsingular solutions.
    complex(wp), allocatable :: solutions(:, :)
      !! The distinct finite nonsingular solutions, shape (n, nsolutions).
    integer :: nreal = 0
      !! Number of real solutions among them.
    real(wp), allocatable :: real_solutions(:, :)
      !! Real parts of the real solutions, shape (n, nreal).
    integer :: nsingular = 0
      !! Number of successful paths ending at finite singular points.
    integer :: ninfinite = 0
      !! Number of successful paths ending at infinity.
    integer :: nfailed = 0
      !! Number of paths the tracker could not follow to t = 0.
    integer :: nduplicates = 0
      !! Number of nonsingular paths ending at an already found
      !! nonsingular solution (indicates path jumping; should be 0).
  end type solve_result

contains

  subroutine solve(sys, res, status, options)
    !! Computes all isolated solutions of the square system `sys` with the
    !! total-degree homotopy. `status` is `err_invalid_input` if the
    !! system is not square or has an equation of degree < 1; failures of
    !! individual paths are reported in `res` instead.

    class(poly_system), intent(inout), target :: sys
    type(solve_result), intent(out) :: res
    integer, intent(out) :: status
    type(solve_options), intent(in), optional :: options

    type(solve_options) :: opts
    type(total_degree_homotopy) :: hom
    type(rng) :: gen
    type(path_info) :: info
    complex(wp), allocatable :: z(:)
    complex(wp), allocatable :: sols(:, :)
    integer :: n
    integer :: k

    if (present(options)) opts = options

    call gen%seed(opts%seed)
    call hom%init(sys, gen, status)
    if (status /= status_ok) return

    n = hom%n
    res%n = n
    res%npaths = hom%nstart
    allocate(res%paths(res%npaths))
    allocate(z(n + 1))
    allocate(sols(n, res%npaths))

    do k = 1, res%npaths
      call hom%start_solution(k, z)
      call track_path(hom, z, 1.0_wp, 0.0_wp, opts%tracker, info)
      associate(p => res%paths(k))
        p%status = info%status
        p%t = info%t
        p%naccepted = info%naccepted
        p%nrejected = info%nrejected
        if (info%status == path_success) then
          call refine(hom, z, opts%refine_maxit)
          call classify(hom, z, opts, p)
          if (p%at_infinity) then
            res%ninfinite = res%ninfinite + 1
          else if (p%singular) then
            res%nsingular = res%nsingular + 1
          else
            call add_solution(p, sols, res, opts%dedup_tol)
          end if
        else
          res%nfailed = res%nfailed + 1
          p%z = z
        end if
      end associate
    end do

    res%solutions = sols(:, 1:res%nsolutions)
    call collect_real(res, opts%real_tol)

  end subroutine solve

  subroutine refine(hom, z, maxit)
    !! Newton refinement of the endpoint on H(z, 0) = 0 (target system and
    !! patch). Stops when the update is at the level of the rounding error
    !! or no longer decreases (e.g. at singular endpoints).

    type(total_degree_homotopy), intent(inout) :: hom
    complex(wp), intent(inout) :: z(:)
    integer, intent(in) :: maxit

    complex(wp) :: h(size(z))
    complex(wp) :: hz(size(z), size(z))
    complex(wp) :: ht(size(z))
    integer :: ipiv(size(z))
    real(wp) :: dnorm
    real(wp) :: dnorm_prev
    integer :: status
    integer :: it

    dnorm_prev = huge(1.0_wp)
    do it = 1, maxit
      call hom%evaluate(z, 0.0_wp, h, hz, ht)
      call lu_factor(hz, ipiv, status)
      if (status /= status_ok) return
      h = -h
      call lu_solve(hz, ipiv, h)
      dnorm = norm_inf(h)
      if (dnorm > dnorm_prev) return
      z = z + h
      if (dnorm <= 4.0_wp*epsilon(1.0_wp)*norm_inf(z)) return
      dnorm_prev = dnorm
    end do

  end subroutine refine

  subroutine classify(hom, z, opts, p)
    !! Stores the endpoint z in p and classifies it (see module
    !! documentation).

    type(total_degree_homotopy), intent(inout) :: hom
    complex(wp), intent(in) :: z(:)
    type(solve_options), intent(in) :: opts
    type(path_result), intent(inout) :: p

    complex(wp) :: h(size(z))
    complex(wp) :: hz(size(z), size(z))
    complex(wp) :: ht(size(z))
    real(wp) :: znorm
    integer :: n

    n = size(z) - 1
    p%z = z
    znorm = norm_inf(z)

    call hom%evaluate(z, 0.0_wp, h, hz, ht)
    p%cond = cond1(hz)
    p%residual = maxval(abs(h(1:n))/znorm**hom%deg)

    p%at_infinity = abs(z(n + 1)) <= opts%infinity_tol*znorm
    p%singular = p%cond > opts%singular_cond
    if (.not. p%at_infinity) then
      p%x = z(1:n)/z(n + 1)
      p%is_real = maxval(abs(aimag(p%x)))  &
        <= opts%real_tol*max(1.0_wp, norm_inf(p%x))
    end if

  end subroutine classify

  subroutine add_solution(p, sols, res, tol)
    !! Adds the finite nonsingular endpoint of p to the list of distinct
    !! solutions, or counts it as a duplicate.

    type(path_result), intent(inout) :: p
    complex(wp), intent(inout) :: sols(:, :)
    type(solve_result), intent(inout) :: res
    real(wp), intent(in) :: tol

    real(wp) :: scale
    integer :: j

    scale = max(1.0_wp, norm_inf(p%x))
    do j = 1, res%nsolutions
      if (norm_inf(p%x - sols(:, j)) <= tol*scale) then
        p%solution_index = j
        res%nduplicates = res%nduplicates + 1
        return
      end if
    end do
    res%nsolutions = res%nsolutions + 1
    sols(:, res%nsolutions) = p%x
    p%solution_index = res%nsolutions

  end subroutine add_solution

  subroutine collect_real(res, tol)
    !! Collects the real parts of the real solutions.

    type(solve_result), intent(inout) :: res
    real(wp), intent(in) :: tol

    logical :: isreal(res%nsolutions)
    integer :: j
    integer :: m

    do j = 1, res%nsolutions
      isreal(j) = maxval(abs(aimag(res%solutions(:, j))))  &
        <= tol*max(1.0_wp, norm_inf(res%solutions(:, j)))
    end do
    res%nreal = count(isreal)
    allocate(res%real_solutions(res%n, res%nreal))
    m = 0
    do j = 1, res%nsolutions
      if (isreal(j)) then
        m = m + 1
        res%real_solutions(:, m) = real(res%solutions(:, j), wp)
      end if
    end do

  end subroutine collect_real

end module polyhomocont__solve
