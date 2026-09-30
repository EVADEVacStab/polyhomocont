# polyhomocont

A Fortran library for computing all isolated complex solutions of
square systems of polynomial equations F(x) = 0, x ∈ Cⁿ, with homotopy
continuation.

This package was written by
[Claude](https://www.anthropic.com/claude) (Anthropic), including the
source code, the test programs, the reference scripts in `reference/`,
and this README.

## Method

The solver uses the total-degree homotopy with the gamma trick,

    H(z, t) = γ t G(z) + (1 - t) F^h(z),   G_i(z) = z_i^{d_i} - b_i w^{d_i},

tracked from t = 1 (start system G, whose d_1 ⋯ d_n solutions are
known) to t = 0 (target system F), with random complex constants γ and
b_i. The paths are followed in projective coordinates
z = (z_1, …, z_n, w) on a random affine patch a·z = 1, so paths
diverging to infinity stay bounded. The path tracker combines a
fourth-order Runge-Kutta predictor with a Newton corrector and adaptive
step size control; steps are only accepted if the corrector converges
quickly, which protects against path jumping.

Near t = 0, the **Cauchy end game** computes the endpoint of each path:
it tracks loops |t| = r in the complex t-plane until the path closes,
which gives the winding number c, and averages the path over the loop
(Cauchy integral formula) at shrinking radii until the estimates agree
and solve the target system. This resolves singular endpoints (multiple
roots) and endpoints at infinity to high accuracy. The endpoints are
then classified as finite or at infinity, regular or singular (winding
number > 1 or ill-conditioned Jacobian), and real or complex; regular
ones are refined with Newton's method.

Positive-dimensional solution components are not treated specially:
paths ending on them are reported as singular endpoints.

## Usage

Systems are passed as objects extending the abstract type
`poly_system`. There are two ways to provide one.

### Sparse representation

Each equation is given by its coefficients and an exponent matrix whose
k-th column is the exponent vector of the k-th term:

```fortran
use polyhomocont

type(sparse_system) :: sys
type(solve_result) :: res
integer :: status

! x^2 + y^2 - 5 = 0,  x y - 2 = 0
call sys%init(2)
call sys%add_equation([1.0_wp, 1.0_wp, -5.0_wp],  &
  reshape([2, 0,  0, 2,  0, 0], [2, 3]), status)
call sys%add_equation([1.0_wp, -2.0_wp],  &
  reshape([1, 1,  0, 0], [2, 2]), status)

call solve(sys, res, status)
print *, res%nsolutions        ! 4
print *, res%solutions         ! complex(wp), shape (n, nsolutions)
print *, res%real_solutions    ! real(wp), shape (n, nreal)
```

Coefficients may be real or complex.

### User-supplied evaluators

A system can also be given through its own evaluator, e.g. when it is
known in closed form. Extend `poly_system` and implement the number of
variables, the total degree of each equation, and a routine that
evaluates F(x) together with the Jacobian dF/dx:

```fortran
type, extends(poly_system) :: my_system
contains
  procedure :: nvars => my_nvars        ! integer: n
  procedure :: degrees => my_degrees    ! integer, allocatable: d(n)
  procedure :: evaluate => my_evaluate  ! (x, f, jac), complex(wp)
end type
```

The degrees must be the exact total degrees of the equations. The
solver needs the homogenized system, which by default is computed from
the affine evaluator as w^{d_i} F_i(z/w); this loses accuracy close to
infinity (w → 0). Types that can evaluate the homogenized polynomials
directly can override `evaluate_homogeneous` (as `sparse_system` does).
See `app/example.f90` for a complete example.

### Results and options

`solve_result` contains

- the distinct finite regular solutions (`solutions`, and the real ones
  in `real_solutions`),
- the distinct finite singular solutions (`singular_solutions`) with
  their `multiplicities` (the number of paths ending there),
- the number of paths going to infinity (`ninfinite`), of failed paths
  (`nfailed`), and of duplicate endpoints of regular paths
  (`nduplicates`, which indicates path jumping and should be zero).

`res%paths(k)` holds the endpoint of every path together with its
classification, winding number, end game accuracy, condition number,
residual and step counts.

The optional argument `options` of type `solve_options` sets the
random seed (results are reproducible for a fixed seed), the
classification tolerances, and in its components `tracker` and
`endgame` the parameters of the path tracker and the end game (e.g.
`endgame%t_start`, the value of t where the end game begins). The end
game can be switched off with `use_endgame = .false.`.

`status` is nonzero (`err_invalid_input`) if the system is not square
or has an equation of degree < 1.

## Building and testing

```sh
fpm build --profile release
fpm test --profile release
fpm run --profile release    # runs app/example.f90
```

Always pass `--profile debug` or `--profile release` explicitly:
without it, fpm compiles with no flags at all, so none of the flags
defined in `fpm.toml` are applied. The code is tested with gfortran and
Intel ifx (`--compiler ifx`).

The tests compare against solutions that are known analytically or were
computed independently with Gröbner bases (sympy) and high-precision
root finding (mpmath); the scripts are in `reference/` (see
`reference/README.md`). Benchmarks with known solution counts:
Katsura-3/4 (8, 16 solutions), cyclic-5 (70 of 120 paths) and cyclic-6
(156 of 720 paths, the others go to infinity).

## Literature

- A. P. Morgan, A. J. Sommese, C. W. Wampler, *Computing singular
  solutions to nonlinear analytic systems*, Numer. Math. 58 (1991) 669
  (Cauchy end game).
- T. L. Lee, T. Y. Li, C. H. Tsai, *HOM4PS-2.0: a software package for
  solving polynomial systems by the polyhedral homotopy continuation
  method*, Computing 83 (2008) 109.
- T. Duff, C. Hill, A. Jensen, K. Lee, A. Leykin, J. Sommars, *Solving
  polynomial systems via homotopy continuation and monodromy*,
  IMA J. Numer. Anal. 39 (2019) 1421,
  [arXiv:1609.08722](https://arxiv.org/abs/1609.08722).

## License

MIT, see [LICENSE](LICENSE).
