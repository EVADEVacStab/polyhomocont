program test_random
!! Tests of the seeded pseudo-random number generator: range, moments,
!! reproducibility and seed dependence.

  use polyhomocont, only : wp
  use polyhomocont__random, only : rng
  use testing, only : assert_true
  use testing, only : assert_small
  use testing, only : finish

  implicit none

  integer, parameter :: nsamples = 200000
  type(rng) :: gen
  type(rng) :: gen2
  real(wp) :: u
  real(wp) :: umin
  real(wp) :: umax
  real(wp) :: mean
  real(wp) :: var
  real(wp) :: absz2
  real(wp) :: dev
  complex(wp) :: z
  complex(wp) :: zmean
  logical :: same
  logical :: differ
  integer :: i

  ! Range and first two moments of the uniform deviates
  call gen%seed(42)
  umin = 1.0_wp
  umax = 0.0_wp
  mean = 0.0_wp
  var = 0.0_wp
  do i = 1, nsamples
    u = gen%uniform()
    umin = min(umin, u)
    umax = max(umax, u)
    mean = mean + u
    var = var + u**2
  end do
  mean = mean/nsamples
  var = var/nsamples - mean**2
  call assert_true("uniform: values in (0, 1)",  &
    umin > 0.0_wp .and. umax < 1.0_wp)
  call assert_small("uniform: mean - 1/2", mean - 0.5_wp, 5.0e-3_wp)
  call assert_small("uniform: variance - 1/12", var - 1.0_wp/12.0_wp,  &
    5.0e-3_wp)

  ! Reproducibility and seed dependence
  call gen%seed(7)
  call gen2%seed(7)
  same = .true.
  do i = 1, 100
    same = same .and. (gen%uniform() == gen2%uniform())
  end do
  call assert_true("seed: equal seeds give equal sequences", same)
  call gen%seed(7)
  call gen2%seed(8)
  differ = .false.
  do i = 1, 100
    differ = differ .or. (gen%uniform() /= gen2%uniform())
  end do
  call assert_true("seed: different seeds give different sequences",  &
    differ)

  ! Unit complex numbers
  call gen%seed(1)
  dev = 0.0_wp
  zmean = (0.0_wp, 0.0_wp)
  do i = 1, nsamples
    z = gen%unit_complex()
    dev = max(dev, abs(abs(z) - 1.0_wp))
    zmean = zmean + z
  end do
  call assert_small("unit_complex: | |z| - 1 |", dev, 1.0e-14_wp)
  call assert_small("unit_complex: |mean|", abs(zmean)/nsamples,  &
    1.0e-2_wp)

  ! Complex normal deviates: E[z] = 0, E[|z|^2] = 1
  zmean = (0.0_wp, 0.0_wp)
  absz2 = 0.0_wp
  do i = 1, nsamples
    z = gen%normal_complex()
    zmean = zmean + z
    absz2 = absz2 + abs(z)**2
  end do
  call assert_small("normal_complex: |mean|", abs(zmean)/nsamples,  &
    1.0e-2_wp)
  call assert_small("normal_complex: E|z|^2 - 1",  &
    absz2/nsamples - 1.0_wp, 1.0e-2_wp)

  call finish()

end program test_random
