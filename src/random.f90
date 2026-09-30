module polyhomocont__random
!! Seeded pseudo-random number generator used for the random constants of
!! the homotopies (gamma trick, random affine patches).
!!
!! The generator is L'Ecuyer's combined multiplicative linear congruential
!! generator (Commun. ACM 31 (1988) 742), period about 2.3e18. All
!! intermediate products fit into 64-bit integers, so the arithmetic is
!! standard-conforming (no reliance on integer overflow). The state lives
!! in an `rng` object rather than in module variables, so the global state
!! of the intrinsic `random_number` of the calling program is untouched,
!! results are reproducible, and independent generators can be used
!! concurrently.

  use, intrinsic :: iso_fortran_env, only : int64
  use polyhomocont__config, only : wp
  use polyhomocont__config, only : pi

  implicit none

  private

  public :: rng

  integer(int64), parameter :: m1 = 2147483563_int64
  integer(int64), parameter :: m2 = 2147483399_int64
  integer(int64), parameter :: a1 = 40014_int64
  integer(int64), parameter :: a2 = 40692_int64

  type :: rng
    !! Pseudo-random number generator. Call `seed` before first use to get
    !! a sequence other than the default one.
    integer(int64) :: s1 = 1234567_int64
    integer(int64) :: s2 = 7654321_int64
  contains
    procedure :: seed => rng_seed
    procedure :: uniform => rng_uniform
    procedure :: unit_complex => rng_unit_complex
    procedure :: normal_complex => rng_normal_complex
  end type rng

contains

  subroutine rng_seed(self, seed)
    !! Initializes the generator from an integer seed. Equal seeds give
    !! equal sequences.

    class(rng), intent(inout) :: self
    integer, intent(in) :: seed

    integer(int64) :: s
    real(wp) :: u
    integer :: i

    s = abs(int(seed, int64))
    self%s1 = mod(s, m1 - 1_int64) + 1_int64
    self%s2 = mod(s*7919_int64 + 104729_int64, m2 - 1_int64) + 1_int64
    ! Discard the first draws to decorrelate nearby seeds
    do i = 1, 16
      u = self%uniform()
    end do

  end subroutine rng_seed

  function rng_uniform(self) result(u)
    !! Uniform deviate in the open interval (0, 1).

    class(rng), intent(inout) :: self
    real(wp) :: u

    integer(int64) :: z

    self%s1 = mod(a1*self%s1, m1)
    self%s2 = mod(a2*self%s2, m2)
    z = self%s1 - self%s2
    if (z < 1_int64) z = z + m1 - 1_int64
    u = real(z, wp)/real(m1, wp)

  end function rng_uniform

  function rng_unit_complex(self) result(z)
    !! Complex number uniformly distributed on the unit circle.

    class(rng), intent(inout) :: self
    complex(wp) :: z

    real(wp) :: phi

    phi = 2.0_wp*pi*self%uniform()
    z = cmplx(cos(phi), sin(phi), kind=wp)

  end function rng_unit_complex

  function rng_normal_complex(self) result(z)
    !! Standard complex normal deviate (real and imaginary parts
    !! independent normal with variance 1/2), via the Box-Muller
    !! transform. A vector of such deviates has a uniformly distributed
    !! direction.

    class(rng), intent(inout) :: self
    complex(wp) :: z

    real(wp) :: r
    real(wp) :: phi

    r = sqrt(-log(self%uniform()))
    phi = 2.0_wp*pi*self%uniform()
    z = cmplx(r*cos(phi), r*sin(phi), kind=wp)

  end function rng_normal_complex

end module polyhomocont__random
