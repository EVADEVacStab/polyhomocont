program test_endgame
!! Tests of the Cauchy end game on total-degree homotopies with singular
!! endpoints: winding numbers and endpoint accuracy for a double and a
!! triple root, and endpoints at infinity.

  use polyhomocont, only : wp
  use polyhomocont, only : sparse_system
  use polyhomocont__random, only : rng
  use polyhomocont__homotopy, only : total_degree_homotopy
  use polyhomocont__tracker, only : tracker_options
  use polyhomocont__tracker, only : path_info
  use polyhomocont__tracker, only : track_path
  use polyhomocont__tracker, only : path_success
  use polyhomocont__endgame, only : endgame_options
  use polyhomocont__endgame, only : endgame_result
  use polyhomocont__endgame, only : cauchy_endgame
  use polyhomocont__endgame, only : endgame_converged
  use testing, only : assert_small
  use testing, only : assert_equal
  use testing, only : assert_true
  use testing, only : finish

  implicit none

  type(rng) :: gen

  call gen%seed(17)
  call check_univariate("double root", [1.0_wp, -1.0_wp, -1.0_wp, 1.0_wp],  &
    1.0_wp, 2, -1.0_wp)
  call check_univariate("triple root",  &
    [1.0_wp, -1.0_wp, -3.0_wp, 5.0_wp, -2.0_wp], 1.0_wp, 3, -2.0_wp)
  call check_infinity()
  call finish()

contains

  subroutine run_paths(sys, z, eg)
    !! Tracks all paths of the total-degree homotopy of sys to t_start and
    !! runs the end game; z(:, k) and eg(k) are the results for path k.

    type(sparse_system), intent(inout), target :: sys
    complex(wp), allocatable, intent(out) :: z(:, :)
    type(endgame_result), allocatable, intent(out) :: eg(:)

    type(total_degree_homotopy) :: hom
    type(tracker_options) :: topts
    type(endgame_options) :: eopts
    type(path_info) :: info
    complex(wp), allocatable :: zk(:)
    integer :: status
    integer :: k

    call hom%init(sys, gen, status)
    allocate(z(sys%nvars() + 1, hom%nstart))
    allocate(eg(hom%nstart))
    allocate(zk(sys%nvars() + 1))
    do k = 1, hom%nstart
      call hom%start_solution(k, zk)
      call track_path(hom, zk, 1.0_wp, eopts%t_start, topts, info)
      call assert_true("tracked to end game", info%status == path_success)
      call cauchy_endgame(hom, zk, eopts, topts, eg(k))
      z(:, k) = eg(k)%z
    end do

  end subroutine run_paths

  subroutine check_univariate(name, coefs, root, mult, other)
    !! Polynomial (x - root)^mult (x - other) with coefficients in
    !! decreasing order: (x - 1)^2 (x + 1) = x^3 - x^2 - x + 1 and
    !! (x - 1)^3 (x + 2) = x^4 - x^3 - 3x^2 + 5x - 2.

    character(len=*), intent(in) :: name
    real(wp), intent(in) :: coefs(:)
    real(wp), intent(in) :: root
    integer, intent(in) :: mult
    real(wp), intent(in) :: other

    type(sparse_system), target :: sys
    complex(wp), allocatable :: z(:, :)
    type(endgame_result), allocatable :: eg(:)
    integer :: exps(1, size(coefs))
    integer :: status
    integer :: deg
    integer :: k
    integer :: nmult
    real(wp) :: err_mult
    real(wp) :: err_other
    logical :: winding_ok
    logical :: converged

    deg = size(coefs) - 1
    exps(1, :) = [(deg - k, k = 0, deg)]
    call sys%init(1)
    call sys%add_equation(coefs, exps, status)
    call run_paths(sys, z, eg)

    nmult = 0
    err_mult = 0.0_wp
    err_other = 0.0_wp
    winding_ok = .true.
    converged = .true.
    do k = 1, size(eg)
      converged = converged .and. eg(k)%status == endgame_converged
      if (abs(z(1, k)/z(2, k) - root) < 0.1_wp) then
        nmult = nmult + 1
        err_mult = max(err_mult, abs(z(1, k)/z(2, k) - root))
        winding_ok = winding_ok .and. eg(k)%winding_number == mult
      else
        err_other = max(err_other, abs(z(1, k)/z(2, k) - other))
        winding_ok = winding_ok .and. eg(k)%winding_number == 1
      end if
    end do
    call assert_equal(name//": paths to the multiple root", nmult, mult)
    call assert_true(name//": winding numbers", winding_ok)
    call assert_true(name//": all end games converged", converged)
    call assert_small(name//": error of the multiple root", err_mult,  &
      1.0e-9_wp)
    call assert_small(name//": error of the simple root", err_other,  &
      1.0e-12_wp)

  end subroutine check_univariate

  subroutine check_infinity()
    !! x y = 1, x^2 = 4: two paths end at the singular point (0 : 1 : 0)
    !! at infinity, i.e. with w = z(3) = 0.

    type(sparse_system), target :: sys
    complex(wp), allocatable :: z(:, :)
    type(endgame_result), allocatable :: eg(:)
    integer :: status
    integer :: k
    integer :: ninf
    real(wp) :: wmax

    call sys%init(2)
    call sys%add_equation([1.0_wp, -1.0_wp],  &
      reshape([1, 1,  0, 0], [2, 2]), status)
    call sys%add_equation([1.0_wp, -4.0_wp],  &
      reshape([2, 0,  0, 0], [2, 2]), status)
    call run_paths(sys, z, eg)

    ninf = 0
    wmax = 0.0_wp
    do k = 1, size(eg)
      if (abs(z(3, k)) < 1.0e-3_wp*maxval(abs(z(:, k)))) then
        ninf = ninf + 1
        wmax = max(wmax, abs(z(3, k))/maxval(abs(z(:, k))))
        ! The point at infinity is (0 : 1 : 0)
        wmax = max(wmax, abs(z(1, k))/maxval(abs(z(:, k))))
      end if
    end do
    call assert_equal("infinity: paths to infinity", ninf, 2)
    call assert_small("infinity: |w|, |z_1| relative to |z|", wmax,  &
      1.0e-9_wp)

  end subroutine check_infinity

end program test_endgame
