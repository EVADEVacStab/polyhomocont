program test_tracker
!! Tests of the path tracker on total-degree homotopies with known
!! solutions: a univariate cubic and a 2x2 system. Checks that all paths
!! arrive at t = 0, that the endpoints are exactly the known solutions,
!! that tracking in two pieces agrees with one piece, and that the step
!! limit is reported.

  use polyhomocont, only : wp
  use polyhomocont, only : sparse_system
  use polyhomocont__random, only : rng
  use polyhomocont__linalg, only : norm_inf
  use polyhomocont__homotopy, only : total_degree_homotopy
  use polyhomocont__tracker, only : tracker_options
  use polyhomocont__tracker, only : path_info
  use polyhomocont__tracker, only : track_path
  use polyhomocont__tracker, only : path_success
  use polyhomocont__tracker, only : path_failed_max_steps
  use testing, only : assert_status
  use testing, only : assert_small
  use testing, only : assert_true
  use testing, only : finish

  implicit none

  type(rng) :: gen

  call gen%seed(3)
  call check_cubic()
  call check_2x2()
  call check_max_steps()
  call finish()

contains

  subroutine check_cubic()
    !! x^3 - 2x + 1 = (x - 1)(x^2 + x - 1): roots 1, (-1 +- sqrt(5))/2.

    type(sparse_system), target :: sys
    type(total_degree_homotopy) :: hom
    type(tracker_options) :: opts
    type(path_info) :: info
    complex(wp) :: z(2)
    real(wp) :: roots(3)
    real(wp) :: err
    logical :: found(3)
    integer :: status
    integer :: k
    integer :: r

    call sys%init(1)
    call sys%add_equation([1.0_wp, -2.0_wp, 1.0_wp],  &
      reshape([3, 1, 0], [1, 3]), status)
    call hom%init(sys, gen, status)
    call assert_status("cubic: init", status)
    roots = [1.0_wp, (-1.0_wp + sqrt(5.0_wp))/2.0_wp,  &
      (-1.0_wp - sqrt(5.0_wp))/2.0_wp]

    found = .false.
    err = 0.0_wp
    do k = 1, hom%nstart
      call hom%start_solution(k, z)
      call track_path(hom, z, 1.0_wp, 0.0_wp, opts, info)
      call assert_true("cubic: path tracked", info%status == path_success)
      r = minloc(abs(z(1)/z(2) - roots), 1)
      found(r) = .true.
      err = max(err, abs(z(1)/z(2) - roots(r)))
    end do
    call assert_true("cubic: all roots found", all(found))
    call assert_small("cubic: max error", err, 1.0e-10_wp)

  end subroutine check_cubic

  subroutine check_2x2()
    !! x^2 + y^2 = 5, x y = 2: solutions (1, 2), (2, 1), (-1, -2),
    !! (-2, -1); the Bezout number 4 equals the number of solutions.

    type(sparse_system), target :: sys
    type(total_degree_homotopy) :: hom
    type(tracker_options) :: opts
    type(path_info) :: info
    type(path_info) :: info2
    complex(wp) :: z(3)
    complex(wp) :: z2(3)
    real(wp) :: sols(2, 4)
    real(wp) :: dist(4)
    real(wp) :: err
    real(wp) :: err_split
    logical :: found(4)
    integer :: status
    integer :: k
    integer :: r

    call sys%init(2)
    call sys%add_equation([1.0_wp, 1.0_wp, -5.0_wp],  &
      reshape([2, 0,  0, 2,  0, 0], [2, 3]), status)
    call sys%add_equation([1.0_wp, -2.0_wp],  &
      reshape([1, 1,  0, 0], [2, 2]), status)
    call hom%init(sys, gen, status)
    call assert_status("2x2: init", status)
    sols = reshape([1.0_wp, 2.0_wp,  2.0_wp, 1.0_wp,  -1.0_wp, -2.0_wp,  &
      -2.0_wp, -1.0_wp], [2, 4])

    found = .false.
    err = 0.0_wp
    err_split = 0.0_wp
    do k = 1, hom%nstart
      call hom%start_solution(k, z)
      z2 = z
      call track_path(hom, z, 1.0_wp, 0.0_wp, opts, info)
      call assert_true("2x2: path tracked", info%status == path_success)
      do r = 1, 4
        dist(r) = norm_inf(z(1:2)/z(3) - sols(:, r))
      end do
      r = minloc(dist, 1)
      found(r) = .true.
      err = max(err, dist(r))

      ! The same path in two pieces
      call track_path(hom, z2, 1.0_wp, 0.5_wp, opts, info2)
      call track_path(hom, z2, 0.5_wp, 0.0_wp, opts, info2)
      err_split = max(err_split, norm_inf(z2 - z)/norm_inf(z))
    end do
    call assert_true("2x2: all solutions found", all(found))
    call assert_small("2x2: max error", err, 1.0e-10_wp)
    call assert_small("2x2: tracking in two pieces", err_split, 1.0e-9_wp)
    call assert_true("2x2: accepted steps counted", info%naccepted > 0)

  end subroutine check_2x2

  subroutine check_max_steps()
    !! With a tiny maximum step size, the step limit must be reported and
    !! t must stay inside (0, 1).

    type(sparse_system), target :: sys
    type(total_degree_homotopy) :: hom
    type(tracker_options) :: opts
    type(path_info) :: info
    complex(wp) :: z(2)
    integer :: status

    call sys%init(1)
    call sys%add_equation([1.0_wp, -2.0_wp],  &
      reshape([2, 0], [1, 2]), status)
    call hom%init(sys, gen, status)
    opts%step_max = 1.0e-3_wp
    opts%step_init = 1.0e-3_wp
    opts%max_steps = 10
    call hom%start_solution(1, z)
    call track_path(hom, z, 1.0_wp, 0.0_wp, opts, info)
    call assert_true("max steps: reported",  &
      info%status == path_failed_max_steps)
    call assert_true("max steps: t in (0, 1)",  &
      info%t > 0.0_wp .and. info%t < 1.0_wp)

  end subroutine check_max_steps

end program test_tracker
