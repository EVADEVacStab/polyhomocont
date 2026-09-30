program test_cyclic
!! Cyclic-n roots problem,
!!   sum_{i=1}^{n} prod_{j=0}^{k-1} x_{i+j} = 0,  k = 1, ..., n - 1,
!!   x_1 x_2 ... x_n - 1 = 0          (indices modulo n),
!! for n = 5 (70 regular solutions) and n = 6 (156 regular solutions),
!! Bjorck and Froberg, J. Symb. Comput. 12 (1991) 329. The Bezout
!! numbers are 120 and 720, so most paths go to infinity; resolving them
!! requires the end game. Besides the counts, the tests check the
!! symmetries of the solution set: invariance under cyclic shifts of the
!! variables and under multiplication by n-th roots of unity.

  use polyhomocont, only : wp
  use polyhomocont, only : sparse_system
  use polyhomocont, only : solve
  use polyhomocont, only : solve_result
  use polyhomocont__config, only : pi
  use testing, only : assert_status
  use testing, only : assert_small
  use testing, only : assert_equal
  use testing, only : set_distance
  use testing, only : finish

  implicit none

  call check_cyclic(5, 70)
  call check_cyclic(6, 156)
  call finish()

contains

  subroutine build_cyclic(n, sys)

    integer, intent(in) :: n
    type(sparse_system), intent(out) :: sys

    integer :: exps(n, n)
    integer :: status
    integer :: k
    integer :: i
    integer :: j

    call sys%init(n)
    do k = 1, n - 1
      exps = 0
      do i = 1, n
        do j = 0, k - 1
          exps(mod(i - 1 + j, n) + 1, i) = 1
        end do
      end do
      call sys%add_equation([(1.0_wp, i = 1, n)], exps, status)
      call assert_status("cyclic: add equation", status)
    end do
    call sys%add_equation([1.0_wp, -1.0_wp],  &
      reshape([[(1, i = 1, n)], [(0, i = 1, n)]], [n, 2]), status)
    call assert_status("cyclic: add equation", status)

  end subroutine build_cyclic

  subroutine check_cyclic(n, nsols)

    integer, intent(in) :: n
    integer, intent(in) :: nsols

    type(sparse_system) :: sys
    type(solve_result) :: res
    complex(wp) :: omega
    character(len=12) :: name
    integer :: status
    integer :: k

    write(name, "(a,i0)") "cyclic-", n
    call build_cyclic(n, sys)
    call solve(sys, res, status)
    call assert_status(trim(name)//": solve", status)
    call assert_equal(trim(name)//": solutions", res%nsolutions, nsols)
    call assert_equal(trim(name)//": singular solutions", res%nsingular, 0)
    call assert_equal(trim(name)//": paths at infinity", res%ninfinite,  &
      res%npaths - nsols)
    call assert_equal(trim(name)//": failed", res%nfailed, 0)
    call assert_equal(trim(name)//": duplicates", res%nduplicates, 0)
    call assert_small(trim(name)//": max residual",  &
      maxval([(res%paths(k)%residual, k = 1, res%npaths)],  &
      mask=[(res%paths(k)%solution_index > 0, k = 1, res%npaths)]),  &
      1.0e-13_wp)

    ! Symmetries of the solution set
    call assert_small(trim(name)//": invariant under cyclic shift",  &
      set_distance(res%solutions, cshift(res%solutions, 1, dim=1)),  &
      1.0e-10_wp)
    omega = cmplx(cos(2.0_wp*pi/n), sin(2.0_wp*pi/n), kind=wp)
    call assert_small(trim(name)//": invariant under x -> omega x",  &
      set_distance(res%solutions, omega*res%solutions), 1.0e-10_wp)

  end subroutine check_cyclic

end program test_cyclic
