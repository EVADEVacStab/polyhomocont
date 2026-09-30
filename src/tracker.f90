module polyhomocont__tracker
!! Predictor-corrector path tracker for homotopies H(z, t) = 0, following
!! a solution path z(t) from t = t_start down to t = t_end.
!!
!! Each step of size dt consists of
!!   - a predictor: one classical fourth-order Runge-Kutta step for the
!!     Davidenko equation dz/dt = -(dH/dz)^-1 dH/dt, and
!!   - a corrector: Newton's method for H(z, t - dt) = 0 at fixed t,
!!     started from the predicted point.
!! A step is accepted if the corrector converges to the relative tolerance
!! `corrector_tol` within `corrector_maxit` iterations and contracts
!! (each Newton update at most `contraction` times the previous one).
!! Requiring fast convergence from a nearby prediction keeps the corrector
!! on the path being followed, which guards against path jumping
!! (HOM4PS-2.0, Lee et al. 2008, Sec. 4). Rejected steps halve the step
!! size; after `grow_after` consecutive accepted steps it is doubled (up
!! to `step_max`). The tracker fails if the step size drops below
!! `step_min` or the number of steps exceeds `max_steps`.

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
  public :: track_path
  public :: newton_correct
  public :: path_success
  public :: path_failed_min_step
  public :: path_failed_max_steps

  ! Path status codes (component `status` of path_info)
  integer, parameter :: path_success = 0
    !! The path was tracked to t_end.
  integer, parameter :: path_failed_min_step = 1
    !! The step size dropped below step_min.
  integer, parameter :: path_failed_max_steps = 2
    !! The maximum number of steps was reached.

  type :: tracker_options
    !! Parameters of the path tracker.
    real(wp) :: step_init = 0.05_wp
      !! Initial step size in t.
    real(wp) :: step_max = 0.1_wp
      !! Largest step size in t.
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
      !! Maximum number of steps (accepted + rejected) per call.
  end type tracker_options

  type :: path_info
    !! Outcome of track_path.
    integer :: status = path_success
      !! path_success or one of the path_failed_* codes.
    real(wp) :: t = 1.0_wp
      !! Value of t reached (t_end on success).
    integer :: naccepted = 0
      !! Number of accepted steps.
    integer :: nrejected = 0
      !! Number of rejected steps.
    real(wp) :: step = 0.0_wp
      !! Step size at the end of tracking.
  end type path_info

contains

  subroutine track_path(hom, z, t_start, t_end, opts, info)
    !! Tracks the solution z of H(z, t_start) = 0 to t = t_end < t_start.
    !! On exit, z is the solution at info%t (= t_end on success, else the
    !! last accepted point).

    class(homotopy), intent(inout) :: hom
    complex(wp), intent(inout) :: z(:)
    real(wp), intent(in) :: t_start
    real(wp), intent(in) :: t_end
    type(tracker_options), intent(in) :: opts
    type(path_info), intent(out) :: info

    complex(wp) :: zp(size(z))
    real(wp) :: t
    real(wp) :: tn
    real(wp) :: dt
    real(wp) :: step
    integer :: nsuccess
    logical :: ok

    t = t_start
    step = min(opts%step_init, opts%step_max)
    nsuccess = 0
    info%naccepted = 0
    info%nrejected = 0

    do while (t > t_end)
      if (info%naccepted + info%nrejected >= opts%max_steps) then
        info%status = path_failed_max_steps
        info%t = t
        info%step = step
        return
      end if

      ! Do not overshoot t_end; land exactly on it
      dt = min(step, t - t_end)
      tn = t - dt
      if (tn - t_end <= 1.0e-3_wp*step) then
        tn = t_end
        dt = t - t_end
      end if

      call predict_rk4(hom, z, t, dt, zp, ok)
      if (ok) call newton_correct(hom, zp, tn, opts%corrector_tol,  &
        opts%corrector_maxit, opts%contraction, ok)

      if (ok) then
        z = zp
        t = tn
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
          info%t = t
          info%step = step
          return
        end if
      end if
    end do

    info%status = path_success
    info%t = t
    info%step = step

  end subroutine track_path

  subroutine predict_rk4(hom, z, t, dt, zp, ok)
    !! One RK4 step of dz/dt = -(dH/dz)^-1 dH/dt from (z, t) to t - dt.
    !! `ok` is false if one of the Jacobians is singular.

    class(homotopy), intent(inout) :: hom
    complex(wp), intent(in) :: z(:)
    real(wp), intent(in) :: t
    real(wp), intent(in) :: dt
    complex(wp), intent(out) :: zp(:)
    logical, intent(out) :: ok

    complex(wp) :: k1(size(z))
    complex(wp) :: k2(size(z))
    complex(wp) :: k3(size(z))
    complex(wp) :: k4(size(z))

    ! The step goes in the direction of decreasing t: increment -dt
    call tangent(hom, z, t, k1, ok)
    if (.not. ok) return
    call tangent(hom, z - 0.5_wp*dt*k1, t - 0.5_wp*dt, k2, ok)
    if (.not. ok) return
    call tangent(hom, z - 0.5_wp*dt*k2, t - 0.5_wp*dt, k3, ok)
    if (.not. ok) return
    call tangent(hom, z - dt*k3, t - dt, k4, ok)
    if (.not. ok) return
    zp = z - dt/6.0_wp*(k1 + 2.0_wp*k2 + 2.0_wp*k3 + k4)

  end subroutine predict_rk4

  subroutine tangent(hom, z, t, dzdt, ok)
    !! Solves dH/dz * dzdt = -dH/dt at (z, t).

    class(homotopy), intent(inout) :: hom
    complex(wp), intent(in) :: z(:)
    real(wp), intent(in) :: t
    complex(wp), intent(out) :: dzdt(:)
    logical, intent(out) :: ok

    complex(wp) :: h(size(z))
    complex(wp) :: hz(size(z), size(z))
    integer :: ipiv(size(z))
    integer :: status

    call hom%evaluate(z, t, h, hz, dzdt)
    call lu_factor(hz, ipiv, status)
    ok = status == status_ok
    if (.not. ok) return
    dzdt = -dzdt
    call lu_solve(hz, ipiv, dzdt)

  end subroutine tangent

  subroutine newton_correct(hom, z, t, tol, maxit, contraction, ok)
    !! Newton's method for H(z, t) = 0 at fixed t, starting from z. `ok`
    !! is true if an update satisfies ||dz|| <= tol*||z|| within maxit
    !! iterations, with every update at most `contraction` times the
    !! previous one; z then holds the corrected point.

    class(homotopy), intent(inout) :: hom
    complex(wp), intent(inout) :: z(:)
    real(wp), intent(in) :: t
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
