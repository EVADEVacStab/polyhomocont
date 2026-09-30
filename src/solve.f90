module polyhomocont__solve
!! Driver that computes all isolated solutions of a square polynomial
!! system with the total-degree homotopy, and post-processes the path
!! endpoints.
!!
!! For each of the d_1 * ... * d_n start solutions, the path is tracked
!! in projective coordinates from t = 1 to the start of the end game
!! (t = 0.1 by default); the Cauchy end game then computes the endpoint
!! z(0) and the winding number c of the path. Each endpoint is classified:
!!   - at infinity if |w| <= infinity_tol * ||z|| (w = z(n+1)), otherwise
!!     finite with x = z(1:n)/w;
!!   - singular if c > 1 or the condition number of the Jacobian of the
!!     (patched, homogenized) target system exceeds singular_cond;
!!     regular endpoints are refined with Newton's method;
!!   - real if ||Im x|| <= real_tol * max(1, ||x||).
!! The distinct finite regular endpoints are collected in `solutions`.
!! Two regular paths cannot end at the same point for a generic
!! homotopy, so such duplicates indicate path jumping and are counted in
!! `nduplicates`. The finite singular endpoints are clustered into the
!! distinct `singular_solutions`; the number of paths ending at each is
!! its multiplicity (for an isolated solution, the total-degree homotopy
!! has as many paths ending at it as its multiplicity).

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
  use polyhomocont__endgame, only : endgame_options
  use polyhomocont__endgame, only : endgame_result
  use polyhomocont__endgame, only : cauchy_endgame
  use polyhomocont__endgame, only : endgame_failed

  implicit none

  private

  public :: solve_options
  public :: path_result
  public :: solve_result
  public :: solve
  public :: path_endgame_failed

  integer, parameter :: path_endgame_failed = 3
    !! Path status: the end game could not compute an endpoint estimate.
    !! (Codes 0-2 are the tracker codes path_success, path_failed_*.)

  type :: solve_options
    !! Parameters of the solver.
    type(tracker_options) :: tracker
      !! Parameters of the path tracker.
    type(endgame_options) :: endgame
      !! Parameters of the Cauchy end game.
    logical :: use_endgame = .true.
      !! If false, paths are tracked directly to t = 0 (singular endpoints
      !! are then not resolved).
    integer :: seed = 1
      !! Seed for the random constants (gamma, start system, patch).
      !! Results are reproducible for a fixed seed.
    integer :: refine_maxit = 5
      !! Maximum number of Newton iterations refining regular endpoints.
    real(wp) :: infinity_tol = 1.0e-8_wp
      !! Endpoints with |w| <= infinity_tol * ||z|| are at infinity.
    real(wp) :: singular_cond = 1.0e10_wp
      !! Endpoints whose Jacobian has a larger condition number are
      !! singular.
    real(wp) :: real_tol = 1.0e-8_wp
      !! Solutions with ||Im x|| <= real_tol * max(1, ||x||) are real.
    real(wp) :: dedup_tol = 1.0e-8_wp
      !! Relative distance below which two regular solutions are
      !! identical.
    real(wp) :: singular_dedup_tol = 1.0e-6_wp
      !! Relative distance below which two singular endpoints belong to
      !! the same solution (singular endpoints are less accurate).
  end type solve_options

  type :: path_result
    !! Endpoint and classification of one path.
    integer :: status = path_success
      !! path_success, a tracker failure code (path_failed_*) for the
      !! path up to the end game, or path_endgame_failed.
    real(wp) :: t = 1.0_wp
      !! Value of t reached on the real axis (0 on success).
    complex(wp), allocatable :: z(:)
      !! Endpoint in projective coordinates (z_1, ..., z_n, w).
    complex(wp), allocatable :: x(:)
      !! Affine endpoint z(1:n)/w; allocated only if the endpoint is
      !! finite.
    logical :: at_infinity = .false.
    logical :: singular = .false.
    logical :: is_real = .false.
    integer :: winding_number = 0
      !! Winding number from the end game (0 without end game).
    integer :: endgame_status = 0
      !! endgame_converged or endgame_not_converged (see
      !! polyhomocont__endgame).
    real(wp) :: endgame_error = 0.0_wp
      !! Relative difference of the last two end game estimates.
    real(wp) :: cond = 0.0_wp
      !! 1-norm condition number of the Jacobian at the endpoint.
    real(wp) :: residual = 0.0_wp
      !! max_i |F^h_i(z)| / ||z||^{d_i}.
    integer :: solution_index = 0
      !! Index into solve_result%solutions (regular finite endpoints).
    integer :: singular_index = 0
      !! Index into solve_result%singular_solutions (singular finite
      !! endpoints).
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
      !! Number of distinct finite regular solutions.
    complex(wp), allocatable :: solutions(:, :)
      !! The distinct finite regular solutions, shape (n, nsolutions).
    integer :: nreal = 0
      !! Number of real solutions among them.
    real(wp), allocatable :: real_solutions(:, :)
      !! Real parts of the real solutions, shape (n, nreal).
    integer :: nsingular = 0
      !! Number of distinct finite singular solutions.
    complex(wp), allocatable :: singular_solutions(:, :)
      !! The distinct finite singular solutions, shape (n, nsingular).
    integer, allocatable :: multiplicities(:)
      !! Number of paths ending at each singular solution.
    integer :: ninfinite = 0
      !! Number of paths ending at infinity.
    integer :: nfailed = 0
      !! Number of paths without an endpoint (tracking or end game
      !! failed).
    integer :: nduplicates = 0
      !! Number of regular paths ending at an already found regular
      !! solution (indicates path jumping; should be 0).
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
    complex(wp), allocatable :: z(:)
    complex(wp), allocatable :: sols(:, :)
    complex(wp), allocatable :: ssols(:, :)
    integer, allocatable :: mult(:)
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
    allocate(ssols(n, res%npaths))
    allocate(mult(res%npaths))

    do k = 1, res%npaths
      call hom%start_solution(k, z)
      call endpoint(hom, z, opts, res%paths(k))
      associate(p => res%paths(k))
        if (p%status /= path_success) then
          res%nfailed = res%nfailed + 1
          p%z = z
        else
          call classify(hom, z, opts, p)
          if (p%at_infinity) then
            res%ninfinite = res%ninfinite + 1
          else if (p%singular) then
            call add_singular(p, ssols, mult, res, opts%singular_dedup_tol)
          else
            call add_solution(p, sols, res, opts%dedup_tol)
          end if
        end if
      end associate
    end do

    res%solutions = sols(:, 1:res%nsolutions)
    res%singular_solutions = ssols(:, 1:res%nsingular)
    res%multiplicities = mult(1:res%nsingular)
    call collect_real(res, opts%real_tol)

  end subroutine solve

  subroutine endpoint(hom, z, opts, p)
    !! Tracks the path starting at z and computes its endpoint at t = 0,
    !! with or without end game; on success, z holds the endpoint.

    type(total_degree_homotopy), intent(inout) :: hom
    complex(wp), intent(inout) :: z(:)
    type(solve_options), intent(in) :: opts
    type(path_result), intent(inout) :: p

    type(path_info) :: info
    type(endgame_result) :: eg
    real(wp) :: t_end

    t_end = 0.0_wp
    if (opts%use_endgame) t_end = opts%endgame%t_start
    call track_path(hom, z, 1.0_wp, t_end, opts%tracker, info)
    p%status = info%status
    p%t = real(info%t, wp)
    p%naccepted = info%naccepted
    p%nrejected = info%nrejected
    if (info%status /= path_success .or. .not. opts%use_endgame) return

    call cauchy_endgame(hom, z, opts%endgame, opts%tracker, eg)
    p%naccepted = p%naccepted + eg%naccepted
    p%nrejected = p%nrejected + eg%nrejected
    if (eg%status == endgame_failed) then
      p%status = path_endgame_failed
      return
    end if
    z = eg%z
    p%t = 0.0_wp
    p%winding_number = eg%winding_number
    p%endgame_status = eg%status
    p%endgame_error = eg%error

  end subroutine endpoint

  subroutine refine(hom, z, maxit)
    !! Newton refinement of the endpoint on H(z, 0) = 0 (target system and
    !! patch). Stops when the update is at the level of the rounding error
    !! or no longer decreases.

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
      call hom%evaluate(z, (0.0_wp, 0.0_wp), h, hz, ht)
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
    !! Refines regular endpoints, stores the endpoint z in p and
    !! classifies it (see module documentation).

    type(total_degree_homotopy), intent(inout) :: hom
    complex(wp), intent(inout) :: z(:)
    type(solve_options), intent(in) :: opts
    type(path_result), intent(inout) :: p

    complex(wp) :: h(size(z))
    complex(wp) :: hz(size(z), size(z))
    complex(wp) :: ht(size(z))
    real(wp) :: znorm
    integer :: n

    n = size(z) - 1

    call hom%evaluate(z, (0.0_wp, 0.0_wp), h, hz, ht)
    p%cond = cond1(hz)
    p%singular = p%winding_number > 1 .or. p%cond > opts%singular_cond
    if (.not. p%singular) then
      call refine(hom, z, opts%refine_maxit)
      call hom%evaluate(z, (0.0_wp, 0.0_wp), h, hz, ht)
      p%cond = cond1(hz)
    end if

    p%z = z
    znorm = norm_inf(z)
    p%residual = maxval(abs(h(1:n))/znorm**hom%deg)
    p%at_infinity = abs(z(n + 1)) <= opts%infinity_tol*znorm
    if (.not. p%at_infinity) then
      p%x = z(1:n)/z(n + 1)
      p%is_real = maxval(abs(aimag(p%x)))  &
        <= opts%real_tol*max(1.0_wp, norm_inf(p%x))
    end if

  end subroutine classify

  subroutine add_solution(p, sols, res, tol)
    !! Adds the finite regular endpoint of p to the list of distinct
    !! solutions, or counts it as a duplicate.

    type(path_result), intent(inout) :: p
    complex(wp), intent(inout) :: sols(:, :)
    type(solve_result), intent(inout) :: res
    real(wp), intent(in) :: tol

    integer :: j

    j = find(p%x, sols(:, 1:res%nsolutions), tol)
    if (j > 0) then
      res%nduplicates = res%nduplicates + 1
    else
      res%nsolutions = res%nsolutions + 1
      sols(:, res%nsolutions) = p%x
      j = res%nsolutions
    end if
    p%solution_index = j

  end subroutine add_solution

  subroutine add_singular(p, ssols, mult, res, tol)
    !! Adds the finite singular endpoint of p to the cluster of an equal
    !! singular solution (increasing its multiplicity) or starts a new
    !! one.

    type(path_result), intent(inout) :: p
    complex(wp), intent(inout) :: ssols(:, :)
    integer, intent(inout) :: mult(:)
    type(solve_result), intent(inout) :: res
    real(wp), intent(in) :: tol

    integer :: j

    j = find(p%x, ssols(:, 1:res%nsingular), tol)
    if (j > 0) then
      mult(j) = mult(j) + 1
    else
      res%nsingular = res%nsingular + 1
      ssols(:, res%nsingular) = p%x
      mult(res%nsingular) = 1
      j = res%nsingular
    end if
    p%singular_index = j

  end subroutine add_singular

  function find(x, pts, tol) result(idx)
    !! Index of the first column of pts within relative distance tol of x
    !! (0 if none).

    complex(wp), intent(in) :: x(:)
    complex(wp), intent(in) :: pts(:, :)
    real(wp), intent(in) :: tol
    integer :: idx

    real(wp) :: scale
    integer :: j

    scale = max(1.0_wp, norm_inf(x))
    do j = 1, size(pts, 2)
      if (norm_inf(x - pts(:, j)) <= tol*scale) then
        idx = j
        return
      end if
    end do
    idx = 0

  end function find

  subroutine collect_real(res, tol)
    !! Collects the real parts of the real regular solutions.

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
