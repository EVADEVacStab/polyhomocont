module polyhomocont__endgame
!! Cauchy end game for computing the endpoint z(0) of a homotopy path,
!! including singular endpoints and endpoints at infinity.
!!
!! Near t = 0, a path z(t) is analytic in s = t^(1/c), where the winding
!! number c is the number of loops around t = 0 after which the path
!! returns to itself (c = 1 for regular endpoints; for singular endpoints
!! c can be larger, and the paths of a cycle converge to the same point).
!! By the Cauchy integral formula, z(0) is the mean of z over the closed
!! loop |t| = r traversed c times; with m equally spaced samples per loop
!! the trapezoid rule converges exponentially in m c
!! (Morgan, Sommese, Wampler, Numer. Math. 58 (1991) 669; Bertini).
!!
!! The end game starts at t = t_start on the real axis, with the path
!! already tracked there. At each radius r it tracks loops |t| = r
!! (m arcs per loop) until the path closes, which gives c and an estimate
!! of z(0). It then moves inward along the real axis to
!! r * radius_factor and repeats, until two consecutive estimates agree
!! to the relative tolerance `tol`. Loops that cannot be tracked (e.g.
!! because the circle is too large and encloses other branch points) are
!! skipped by moving inward. In projective coordinates on an affine
!! patch, endpoints at infinity (w = 0) are handled exactly like finite
!! ones.

  use polyhomocont__config, only : wp
  use polyhomocont__config, only : pi
  use polyhomocont__linalg, only : norm_inf
  use polyhomocont__homotopy, only : homotopy
  use polyhomocont__tracker, only : tracker_options
  use polyhomocont__tracker, only : path_info
  use polyhomocont__tracker, only : track_segment
  use polyhomocont__tracker, only : line_segment
  use polyhomocont__tracker, only : arc_segment
  use polyhomocont__tracker, only : path_success

  implicit none

  private

  public :: endgame_options
  public :: endgame_result
  public :: cauchy_endgame
  public :: endgame_converged
  public :: endgame_not_converged
  public :: endgame_failed

  ! End game status codes (component `status` of endgame_result)
  integer, parameter :: endgame_converged = 0
    !! Two consecutive estimates agree to the tolerance.
  integer, parameter :: endgame_not_converged = 1
    !! radius_min was reached before convergence; the estimate with the
    !! smallest difference to its predecessor is returned.
  integer, parameter :: endgame_failed = 2
    !! No estimate could be computed.

  type :: endgame_options
    !! Parameters of the Cauchy end game.
    real(wp) :: t_start = 0.1_wp
      !! Value of t at which the end game starts.
    integer :: nsamples = 8
      !! Number of sample points (arcs) per loop.
    real(wp) :: radius_factor = 0.25_wp
      !! Factor by which the radius shrinks between estimates.
    real(wp) :: radius_min = 1.0e-13_wp
      !! Smallest radius; below it, the end game gives up.
    integer :: max_winding = 16
      !! Largest winding number considered.
    real(wp) :: tol = 1.0e-11_wp
      !! Relative agreement of consecutive estimates for convergence.
    real(wp) :: closure_factor = 0.1_wp
      !! A loop is closed if the distance of its end point to its start
      !! point is below closure_factor times the smallest distance
      !! between consecutive samples (scale-free test, see loops).
  end type endgame_options

  type :: endgame_result
    !! Outcome of cauchy_endgame.
    integer :: status = endgame_failed
    complex(wp), allocatable :: z(:)
      !! Estimate of the endpoint z(0).
    integer :: winding_number = 0
    real(wp) :: error = huge(1.0_wp)
      !! Relative difference between the returned estimate and the
      !! previous one (an estimate of its accuracy).
    real(wp) :: radius = 0.0_wp
      !! Radius at which the returned estimate was computed.
    integer :: naccepted = 0
      !! Accepted tracker steps (loops and radial segments).
    integer :: nrejected = 0
      !! Rejected tracker steps.
  end type endgame_result

contains

  subroutine cauchy_endgame(hom, z, opts, topts, res)
    !! Runs the Cauchy end game for the path through z at t = t_start
    !! (on the real axis). On exit, res%z holds the estimate of z(0) and
    !! z the last point tracked on the real axis.

    class(homotopy), intent(inout) :: hom
    complex(wp), intent(inout) :: z(:)
    type(endgame_options), intent(in) :: opts
    type(tracker_options), intent(in) :: topts
    type(endgame_result), intent(out) :: res

    type(tracker_options) :: segopts
    type(path_info) :: info
    complex(wp) :: est(size(z))
    complex(wp) :: prev(size(z))
    real(wp) :: r
    real(wp) :: rnew
    real(wp) :: err
    integer :: c
    logical :: have_prev
    logical :: ok

    ! Arcs and radial segments are short: allow a single step for each,
    ! refined adaptively if the corrector does not converge
    segopts = topts
    segopts%step_init = 1.0_wp
    segopts%step_max = 1.0_wp

    res%status = endgame_failed
    res%z = z
    r = opts%t_start
    have_prev = .false.

    do
      call loops(hom, z, r, opts, segopts, est, c, ok, res)
      if (ok) then
        if (have_prev) then
          err = norm_inf(est - prev)/max(norm_inf(est), tiny(1.0_wp))
          if (err < res%error) then
            res%status = endgame_not_converged
            res%z = est
            res%winding_number = c
            res%error = err
            res%radius = r
          end if
          if (err <= opts%tol) then
            res%status = endgame_converged
            return
          end if
        end if
        prev = est
        have_prev = .true.
      else
        have_prev = .false.
      end if

      rnew = opts%radius_factor*r
      if (rnew < opts%radius_min) return
      call track_segment(hom, z, line_segment(cmplx(r, 0.0_wp, kind=wp),  &
        cmplx(rnew, 0.0_wp, kind=wp)), segopts, info)
      res%naccepted = res%naccepted + info%naccepted
      res%nrejected = res%nrejected + info%nrejected
      if (info%status /= path_success) return
      r = rnew
    end do

  end subroutine cauchy_endgame

  subroutine loops(hom, z, r, opts, segopts, est, c, ok, res)
    !! Tracks loops |t| = r starting and ending at t = r, with nsamples
    !! arcs per loop, until the path closes. Returns the winding number c
    !! and the Cauchy estimate est = mean of the samples; `ok` is false if
    !! an arc fails or the path does not close within max_winding loops.
    !! z is not modified.
    !!
    !! Closure test: after c loops, the path is closed if the distance of
    !! the current point to the start point is at most closure_factor
    !! times the smallest distance between consecutive samples (or at
    !! the level of the corrector tolerance). If the path has not closed,
    !! it is on another sheet of its cycle, whose distance to the start
    !! point is of the same order as the sample spacing or larger.

    class(homotopy), intent(inout) :: hom
    complex(wp), intent(in) :: z(:)
    real(wp), intent(in) :: r
    type(endgame_options), intent(in) :: opts
    type(tracker_options), intent(in) :: segopts
    complex(wp), intent(out) :: est(:)
    integer, intent(out) :: c
    logical, intent(out) :: ok
    type(endgame_result), intent(inout) :: res

    type(path_info) :: info
    complex(wp) :: zc(size(z))
    complex(wp) :: zprev(size(z))
    complex(wp) :: zsum(size(z))
    real(wp) :: dtheta
    real(wp) :: spacing
    real(wp) :: floor
    integer :: nsamp
    integer :: j

    ok = .false.
    dtheta = 2.0_wp*pi/real(opts%nsamples, wp)
    floor = 10.0_wp*segopts%corrector_tol*norm_inf(z)
    zc = z
    zsum = (0.0_wp, 0.0_wp)
    nsamp = 0
    spacing = huge(1.0_wp)

    do c = 1, opts%max_winding
      do j = 1, opts%nsamples
        zsum = zsum + zc
        nsamp = nsamp + 1
        zprev = zc
        call track_segment(hom, zc,  &
          arc_segment(r, (j - 1)*dtheta, dtheta), segopts, info)
        res%naccepted = res%naccepted + info%naccepted
        res%nrejected = res%nrejected + info%nrejected
        if (info%status /= path_success) return
        spacing = min(spacing, norm_inf(zc - zprev))
      end do
      if (norm_inf(zc - z) <= max(opts%closure_factor*spacing, floor))  &
        then
        est = zsum/real(nsamp, wp)
        ok = .true.
        return
      end if
    end do

  end subroutine loops

end module polyhomocont__endgame
