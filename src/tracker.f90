module polyhomocont__tracker
!! Predictor-corrector path tracker for homotopies H(z, t) = 0.
!!
!! A solution z is followed along a segment t(s), s in [0, 1], of the
!! complex t-plane: a straight line from t0 to t1 (`line_segment`) or an
!! arc t = r exp(i (theta0 + s dtheta)) of the circle |t| = r
!! (`arc_segment`, used by the end game). Each step of size ds consists of
!!   - a predictor: one classical fourth-order Runge-Kutta step for the
!!     Davidenko equation dz/ds = -(dH/dz)^-1 dH/dt dt/ds, and
!!   - a corrector: Newton's method for H(z, t(s + ds)) = 0 at fixed t,
!!     started from the predicted point.
!! A step is accepted if the corrector converges to the relative tolerance
!! `corrector_tol` within `corrector_maxit` iterations and contracts
!! (each Newton update at most `contraction` times the previous one).
!! Requiring fast convergence from a nearby prediction keeps the corrector
!! on the path being followed, which guards against path jumping
!! (HOM4PS-2.0, Lee et al. 2008, Sec. 4). Rejected steps halve the step
!! size; after `grow_after` consecutive accepted steps it is doubled (up
!! to `step_max`). The tracker fails if the step size drops below
!! `step_min` or the number of steps exceeds `max_steps`. Step sizes are
!! fractions of the segment (s units).

  use polyhomocont__config, only : wp
  use polyhomocont__config, only : status_ok
  use polyhomocont__linalg, only : lu_factor
  use polyhomocont__linalg, only : lu_solve
  use polyhomocont__linalg, only : norm_inf
  use polyhomocont__homotopy, only : homotopy

  implicit none

  private

  public :: tracker_options
  public :: path_info
  public :: segment
  public :: line_segment
  public :: arc_segment
  public :: track_segment
  public :: track_path
  public :: path_success
  public :: path_failed_min_step
  public :: path_failed_max_steps

  ! Path status codes (component `status` of path_info)
  integer, parameter :: path_success = 0
    !! The path was tracked to the end of the segment.
  integer, parameter :: path_failed_min_step = 1
    !! The step size dropped below step_min.
  integer, parameter :: path_failed_max_steps = 2
    !! The maximum number of steps was reached.

  integer, parameter :: kind_line = 1
  integer, parameter :: kind_arc = 2

  type :: tracker_options
    !! Parameters of the path tracker. Step sizes are fractions of the
    !! tracked segment; for the main path from t = 1 to the start of the
    !! end game they are (up to the factor 1 - t_endgame) steps in t.
    real(wp) :: step_init = 0.05_wp
      !! Initial step size.
    real(wp) :: step_max = 0.1_wp
      !! Largest step size.
    real(wp) :: step_min = 1.0e-14_wp
      !! Smallest step size; below it, tracking fails.
    real(wp) :: corrector_tol = 1.0e-10_wp
      !! Relative tolerance ||dz|| <= corrector_tol * ||z|| (maximum norm)
      !! for the Newton updates of the corrector.
    integer :: corrector_maxit = 3
      !! Maximum number of Newton iterations per step.
    real(wp) :: contraction = 0.5_wp
      !! Required contraction ||dz_k|| <= contraction * ||dz_{k-1}||.
    integer :: grow_after = 5
      !! Number of consecutive accepted steps before the step size grows.
    integer :: max_steps = 20000
      !! Maximum number of steps (accepted + rejected) per segment.
  end type tracker_options

  type :: segment
    !! A segment t(s), s in [0, 1], of the complex t-plane. Construct with
    !! line_segment or arc_segment.
    integer :: kind = kind_line
    complex(wp) :: t0 = (1.0_wp, 0.0_wp)
      !! Line: start point.
    complex(wp) :: t1 = (0.0_wp, 0.0_wp)
      !! Line: end point.
    real(wp) :: radius = 0.0_wp
      !! Arc: radius r.
    real(wp) :: theta0 = 0.0_wp
      !! Arc: start angle.
    real(wp) :: dtheta = 0.0_wp
      !! Arc: angle swept (signed).
  end type segment

  type :: path_info
    !! Outcome of track_segment.
    integer :: status = path_success
      !! path_success or one of the path_failed_* codes.
    real(wp) :: s = 0.0_wp
      !! Segment parameter reached (1 on success).
    complex(wp) :: t = (0.0_wp, 0.0_wp)
      !! Value of t reached, t(s).
    integer :: naccepted = 0
      !! Number of accepted steps.
    integer :: nrejected = 0
      !! Number of rejected steps.
    real(wp) :: step = 0.0_wp
      !! Step size (s units) proposed for continuing the tracking; can be
      !! used as step_init for a following segment of similar length.
  end type path_info

contains

  pure function line_segment(t0, t1) result(seg)
    !! Straight segment from t0 to t1.

    complex(wp), intent(in) :: t0
    complex(wp), intent(in) :: t1
    type(segment) :: seg

    seg%kind = kind_line
    seg%t0 = t0
    seg%t1 = t1

  end function line_segment

  pure function arc_segment(radius, theta0, dtheta) result(seg)
    !! Arc t = radius * exp(i (theta0 + s dtheta)), s in [0, 1].

    real(wp), intent(in) :: radius
    real(wp), intent(in) :: theta0
    real(wp), intent(in) :: dtheta
    type(segment) :: seg

    seg%kind = kind_arc
    seg%radius = radius
    seg%theta0 = theta0
    seg%dtheta = dtheta

  end function arc_segment

  pure function seg_t(seg, s) result(t)
    !! Point t(s) of the segment. The end points are returned exactly.

    type(segment), intent(in) :: seg
    real(wp), intent(in) :: s
    complex(wp) :: t

    real(wp) :: theta

    if (seg%kind == kind_line) then
      if (s >= 1.0_wp) then
        t = seg%t1
      else
        t = seg%t0 + s*(seg%t1 - seg%t0)
      end if
    else
      theta = seg%theta0 + s*seg%dtheta
      t = seg%radius*cmplx(cos(theta), sin(theta), kind=wp)
    end if

  end function seg_t

  pure function seg_dt(seg, s) result(dt)
    !! Derivative dt/ds of the segment.

    type(segment), intent(in) :: seg
    real(wp), intent(in) :: s
    complex(wp) :: dt

    if (seg%kind == kind_line) then
      dt = seg%t1 - seg%t0
    else
      dt = cmplx(0.0_wp, seg%dtheta, kind=wp)*seg_t(seg, s)
    end if

  end function seg_dt

  subroutine track_path(hom, z, t_start, t_end, opts, info)
    !! Tracks the solution z of H(z, t_start) = 0 along the real axis to
    !! t = t_end.

    class(homotopy), intent(inout) :: hom
    complex(wp), intent(inout) :: z(:)
    real(wp), intent(in) :: t_start
    real(wp), intent(in) :: t_end
    type(tracker_options), intent(in) :: opts
    type(path_info), intent(out) :: info

    call track_segment(hom, z,  &
      line_segment(cmplx(t_start, 0.0_wp, kind=wp),  &
      cmplx(t_end, 0.0_wp, kind=wp)), opts, info)

  end subroutine track_path

  subroutine track_segment(hom, z, seg, opts, info)
    !! Tracks the solution z of H(z, t(0)) = 0 along the segment to t(1).
    !! On exit, z is the solution at info%t (= t(1) on success, else the
    !! last accepted point).

    class(homotopy), intent(inout) :: hom
    complex(wp), intent(inout) :: z(:)
    type(segment), intent(in) :: seg
    type(tracker_options), intent(in) :: opts
    type(path_info), intent(out) :: info

    complex(wp) :: zp(size(z))
    real(wp) :: s
    real(wp) :: sn
    real(wp) :: ds
    real(wp) :: step
    integer :: nsuccess
    logical :: ok

    s = 0.0_wp
    step = min(opts%step_init, opts%step_max)
    nsuccess = 0
    info%naccepted = 0
    info%nrejected = 0
    info%status = path_success

    do while (s < 1.0_wp)
      if (info%naccepted + info%nrejected >= opts%max_steps) then
        info%status = path_failed_max_steps
        exit
      end if

      ! Do not overshoot the end of the segment; land exactly on it
      ds = min(step, 1.0_wp - s)
      sn = s + ds
      if (1.0_wp - sn <= 1.0e-3_wp*step) then
        sn = 1.0_wp
        ds = 1.0_wp - s
      end if

      call predict_rk4(hom, seg, z, s, ds, zp, ok)
      if (ok) call newton_correct(hom, zp, seg_t(seg, sn),  &
        opts%corrector_tol, opts%corrector_maxit, opts%contraction, ok)

      if (ok) then
        z = zp
        s = sn
        info%naccepted = info%naccepted + 1
        nsuccess = nsuccess + 1
        if (nsuccess >= opts%grow_after) then
          step = min(2.0_wp*step, opts%step_max)
          nsuccess = 0
        end if
      else
        info%nrejected = info%nrejected + 1
        nsuccess = 0
        step = 0.5_wp*step
        if (step < opts%step_min) then
          info%status = path_failed_min_step
          exit
        end if
      end if
    end do

    info%s = s
    info%t = seg_t(seg, s)
    info%step = step

  end subroutine track_segment

  subroutine predict_rk4(hom, seg, z, s, ds, zp, ok)
    !! One RK4 step of dz/ds = -(dH/dz)^-1 dH/dt dt/ds from (z, s) to
    !! s + ds. `ok` is false if one of the Jacobians is singular.

    class(homotopy), intent(inout) :: hom
    type(segment), intent(in) :: seg
    complex(wp), intent(in) :: z(:)
    real(wp), intent(in) :: s
    real(wp), intent(in) :: ds
    complex(wp), intent(out) :: zp(:)
    logical, intent(out) :: ok

    complex(wp) :: k1(size(z))
    complex(wp) :: k2(size(z))
    complex(wp) :: k3(size(z))
    complex(wp) :: k4(size(z))

    call tangent(hom, seg, z, s, k1, ok)
    if (.not. ok) return
    call tangent(hom, seg, z + 0.5_wp*ds*k1, s + 0.5_wp*ds, k2, ok)
    if (.not. ok) return
    call tangent(hom, seg, z + 0.5_wp*ds*k2, s + 0.5_wp*ds, k3, ok)
    if (.not. ok) return
    call tangent(hom, seg, z + ds*k3, s + ds, k4, ok)
    if (.not. ok) return
    zp = z + ds/6.0_wp*(k1 + 2.0_wp*k2 + 2.0_wp*k3 + k4)

  end subroutine predict_rk4

  subroutine tangent(hom, seg, z, s, dzds, ok)
    !! Solves dH/dz * dzds = -dH/dt * dt/ds at (z, t(s)).

    class(homotopy), intent(inout) :: hom
    type(segment), intent(in) :: seg
    complex(wp), intent(in) :: z(:)
    real(wp), intent(in) :: s
    complex(wp), intent(out) :: dzds(:)
    logical, intent(out) :: ok

    complex(wp) :: h(size(z))
    complex(wp) :: hz(size(z), size(z))
    integer :: ipiv(size(z))
    integer :: status

    call hom%evaluate(z, seg_t(seg, s), h, hz, dzds)
    call lu_factor(hz, ipiv, status)
    ok = status == status_ok
    if (.not. ok) return
    dzds = -dzds*seg_dt(seg, s)
    call lu_solve(hz, ipiv, dzds)

  end subroutine tangent

  subroutine newton_correct(hom, z, t, tol, maxit, contraction, ok)
    !! Newton's method for H(z, t) = 0 at fixed t, starting from z. `ok`
    !! is true if an update satisfies ||dz|| <= tol*||z|| within maxit
    !! iterations, with every update at most `contraction` times the
    !! previous one; z then holds the corrected point.

    class(homotopy), intent(inout) :: hom
    complex(wp), intent(inout) :: z(:)
    complex(wp), intent(in) :: t
    real(wp), intent(in) :: tol
    integer, intent(in) :: maxit
    real(wp), intent(in) :: contraction
    logical, intent(out) :: ok

    complex(wp) :: h(size(z))
    complex(wp) :: hz(size(z), size(z))
    complex(wp) :: ht(size(z))
    integer :: ipiv(size(z))
    integer :: status
    real(wp) :: dnorm
    real(wp) :: dnorm_prev
    integer :: it

    ok = .false.
    dnorm_prev = huge(1.0_wp)
    do it = 1, maxit
      call hom%evaluate(z, t, h, hz, ht)
      call lu_factor(hz, ipiv, status)
      if (status /= status_ok) return
      h = -h
      call lu_solve(hz, ipiv, h)
      z = z + h
      dnorm = norm_inf(h)
      if (dnorm <= tol*norm_inf(z)) then
        ok = .true.
        return
      end if
      if (it > 1 .and. dnorm > contraction*dnorm_prev) return
      dnorm_prev = dnorm
    end do

  end subroutine newton_correct

end module polyhomocont__tracker
