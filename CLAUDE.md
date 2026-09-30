# CLAUDE.md

Guidance for Claude Code when working on **polyhomocont**.

## Goal

A Fortran library (fpm package) that finds **all isolated complex
solutions** of a square system of polynomial equations
F(x) = 0, x ∈ Cⁿ, using **homotopy continuation**. Real solutions
are then picked out among the complex ones. The library should be
robust (no lost paths, correct handling of solutions at infinity and
singular solutions), be usable from other Fortran projects, and be
tested against independently known solutions.

## Directory layout

```
../                   (/home/thomas/Work/Tools/polyhomocont)
├── polyhomocont/     this fpm project (git repo); all new code goes here
├── fortransitions/   reference project for code style, structure and
│                     fpm/compiler configuration (read-only; do not edit)
└── Notes/            literature (read-only)
    ├── lee2008.pdf        Lee, Li, Tsai: HOM4PS-2.0, polyhedral homotopy
    │                      (mixed cells, polyhedral-linear homotopy,
    │                      curve jumping, scaling, end games)
    └── 1609.08722v4.pdf   Duff et al.: solving polynomial systems via
                           homotopy continuation and monodromy
```

Read PDFs with `pdftotext file.pdf -` (poppler-utils is installed).
Other packages by the same author (e.g. `gradmin`, `odeint`) live in
`/home/thomas/Work/Tools/` and can be used as dependencies via git URLs
the same way `fortransitions/fpm.toml` does.

## Build and test

```sh
fpm build --profile debug
fpm test  --profile debug      # during development (bounds/FPE checks)
fpm test  --profile release
fpm run   --profile release
```

Always pass `--profile debug` or `--profile release`: without a profile
fpm applies none of the flags from `fpm.toml`. The `fpm.toml` should
mirror `fortransitions/fpm.toml` (`[build]`, `[install]`, `[fortran]`
with `implicit-typing = false`, `[features]` with the gfortran/ifx debug
and release flags, `[profiles]`). Compilers: gfortran 15.2 (default) and
ifx 2025.3; code must compile warning-free with both. For ifx keep
`-fp-model=precise`. Build with ifx via
`fpm build --compiler ifx --profile ...`; ifx is only available after
`source ~/intel/oneapi/setvars.sh` (check `echo $SETVARS_COMPLETED`
or `which ifx`; if missing, ask the user to source it). Intel MKL is
installed under `~/intel/oneapi/mkl` (relevant if LAPACK is used).

## Code style (follow fortransitions)

- Free-form Fortran 2008+, lowercase keywords, 2-space indentation,
  lines ≤ 79 characters, continuation as `  &` at the end of the line.
- One module per file in `src/`. Modules are named
  `polyhomocont__<topic>` (double underscore, `module-naming = false`);
  the top-level module `polyhomocont` only re-exports the public API.
- A `polyhomocont__config` module defines `wp = real64`, constants and
  integer status codes (`status_ok = 0`, `err_*`). Errors are reported
  via `integer, intent(out) :: status` arguments, never by stopping.
- Every module: `implicit none`, `private` by default, explicit
  `public ::` statements one per line. One symbol per `use ..., only :`
  line.
- One variable per declaration line; all dummy arguments carry
  `intent`. Functions use `result(...)`. Blank line after the
  declaration block and before `end`. Named `end function foo` etc.
- Optional arguments: copy into a local `name_` with a default, then
  `if (present(name)) name_ = name`.
- Literals with kind suffix (`1.0_wp`, `(0.0_wp, 1.0_wp)`); no implicit
  kind conversions.
- FORD-style `!!` doc comments directly below every module, type,
  component and public procedure; plain `!` for inline comments.
  Explain the numerics (which equation / which paper section) in the
  doc comments.
- User-supplied evaluators (F, dF/dx) via abstract derived types with
  deferred type-bound procedures, as in `cosmotransitions__potentials`.

## Tests

- Each test is a standalone program in `test/` (auto-discovered), one
  per module plus problem-level tests (`solve.f90`, `katsura.f90`).
  Shared helpers live in `test/testing.f90` (module `testing`):
  `assert_status`, `assert_close`, `assert_small`, `assert_true`,
  `assert_equal` print `ok   ...` / `FAIL ...` and count failures;
  `call finish()` ends each test (`error stop 1` on failure). It also
  has finite-difference Jacobian checks, a homogeneity check and
  `set_distance` for comparing solution sets up to ordering. Test
  modules must not use nested `contains` inside internal procedures.
- Before each commit, run the tests in all four configurations
  (gfortran/ifx x debug/release) and check for compiler warnings, e.g.
  `for c in gfortran ifx; do for p in debug release; do
  fpm test --compiler $c --profile $p; done; done` (use
  `fpm clean --all` first to see warnings of cached objects).
- Verification, in order of preference:
  1. **No external tool**: systems with closed-form solutions, and
     benchmark families with known root counts (Katsura-n: 2ⁿ,
     cyclic-5: 70, cyclic-6: 156); systems with solutions at infinity
     and singular solutions. Self-consistency checks: residuals,
     condition numbers, identical solution sets for several random γ.
  2. **Python (primary reference tool)**: `numpy.roots` for univariate
     tests; `sympy` (`groebner`, `solve_poly_system`) for small
     systems (n ≲ 4, low degree). Gröbner bases are an independent
     method, so they share no failure modes with continuation; evaluate
     the exact solutions to high precision for the references. Also
     gives multiplicities.
  3. **HomotopyContinuation.jl (Julia)** for larger systems, solution
     counts, and mixed volumes (via MixedSubdivisions.jl, for testing
     the polyhedral homotopy later). Julia scripts directly, or from
     Python via `juliacall`.
  Not used for now: phcpy (needs building PHCpack with GNAT/Ada),
  Bertini 1.x (file-driven CLI, no Python interface), pybertini
  (immature). Later option: α-theory certification of solutions
  (Duff et al. Sec. 5.3).
- Python: use the pyenv virtualenv `polyhomocont` (Python 3.13.11,
  activated via `PYENV_VERSION=polyhomocont`; check with
  `pyenv version`). Install reference-script dependencies (numpy,
  sympy, mpmath, ...) only into this virtualenv, never system-wide.
  Julia, phcpy and Bertini are not installed; ask the user before
  installing anything outside the virtualenv.
- Keep the generating scripts in `reference/` with a README (tool
  versions, how to run), as in `fortransitions/reference/`.
- Solutions come out in no particular order: store references as lists
  of solutions and match computed and reference solution sets up to
  permutation within a tolerance (also checking that the counts
  agree).

## Algorithmic roadmap (initial plan, revise as work progresses)

1. **System interface and evaluation**: systems are passed as objects
   of an abstract derived type (e.g. `poly_system`) with deferred
   type-bound procedures that evaluate F(x) and dF/dx (ideally one
   combined procedure returning both) and report the dimension n and
   the degree d_i of each equation (needed for the start system and
   Bézout count). Two ways to use it:
   - a concrete built-in extension holding a sparse representation
     (support exponent matrix + complex coefficients per equation),
     with efficient joint evaluation of F and dF/dx (shared power
     tables);
   - user-supplied extensions implementing their own evaluator for F
     and dF/dx (e.g. for systems in closed form, as in
     `cosmotransitions__potentials`).
   The tracker must only work through the abstract interface, never
   assume the sparse representation. Features that need the monomial
   supports (polyhedral homotopy, coefficient scaling) may require the
   sparse type or an optional procedure providing the supports.
2. **Linear algebra**: complex LU with partial pivoting for the Newton
   and predictor steps; condition number estimate for singularity
   detection.
3. **Total-degree (linear) homotopy** with the gamma trick:
   H(x, t) = (1 - t)·γ·G(x) + t·F(x), G_i = x_i^{d_i} - 1, random
   complex γ. Start solutions are roots of unity.
4. **Path tracker**: predictor (Euler / RK4, possibly Hermite) +
   Newton corrector, adaptive step size control, guarding against path
   jumping (HOM4PS-2.0 Sec. 4). Work in projective coordinates or on a
   random affine patch so paths diverging to infinity stay bounded.
5. **End game** near t = 1 (power-series or Cauchy end game) for
   singular endpoints and winding numbers (HOM4PS-2.0 Sec. 5.3).
6. **Post-processing**: Newton refinement, classification
   (finite/infinite, regular/singular, real/complex), deduplication,
   multiplicities.
7. **Later**: coefficient scaling (HOM4PS Sec. 5.2), polyhedral homotopy
   via mixed volume / mixed cells (HOM4PS Sec. 2–3), parameter homotopy
   and monodromy solving (Duff et al.).

## Architecture (first milestone)

| File              | Module                    | Content |
| ----------------- | ------------------------- | ------- |
| `config.f90`      | `polyhomocont__config`    | `wp`, constants, status codes |
| `random.f90`      | `polyhomocont__random`    | own PRNG object (L'Ecuyer combined MLCG), no use of the global `random_number` state |
| `linalg.f90`      | `polyhomocont__linalg`    | complex LU with partial pivoting, solve, 1-norm condition number |
| `system.f90`      | `polyhomocont__system`    | abstract `poly_system` (user extension point) |
| `sparse.f90`      | `polyhomocont__sparse`    | `sparse_system`: built-in sparse polynomial system |
| `homotopy.f90`    | `polyhomocont__homotopy`  | abstract `homotopy` + `total_degree_homotopy` |
| `tracker.f90`     | `polyhomocont__tracker`   | predictor-corrector path tracker, `tracker_options` |
| `solve.f90`       | `polyhomocont__solve`     | `solve` driver, post-processing, result types |
| `polyhomocont.f90`| `polyhomocont`            | re-exports the public API |

Conventions:
- `poly_system` has deferred `nvars`, `degrees` and `evaluate(x, f, jac)`
  (affine F and dF/dx together). `evaluate_homogeneous(z, f, jac)` has a
  default implementation via the affine evaluator,
  F^h_i(z) = w^{d_i} F_i(z(1:n)/w), w = z(n+1), which loses accuracy
  as w -> 0; `sparse_system` overrides it with the exact homogenized
  polynomial.
- Projective coordinates: z(1:n+1), homogenizing coordinate **last**
  (w = z(n+1), x = z(1:n)/w). Tracking happens on a random affine patch
  a·z = 1, giving a square (n+1)x(n+1) system.
- Homotopy parameter runs from t = 1 (start system) to t = 0 (target),
  H(z, t) = γ t G(z) + (1 - t) F^h(z), G_i = z_i^{d_i} - w^{d_i}, so
  that the end game region near t = 0 has full floating-point
  resolution (convention of Bertini / HomotopyContinuation.jl).
- The tracker only sees `class(homotopy)`, so other homotopies
  (polyhedral, parameter) can be added later.
- Randomness (γ, patch) comes from a seeded PRNG object, so results are
  reproducible; the seed is an option.

## Design decisions (settled with the user)

- **Linear algebra**: self-contained complex LU with partial pivoting,
  no LAPACK dependency (systems are small). A LAPACK/MKL backend may be
  added later as an option.
- **Parallelisation**: none for now, but keep it possible: every path
  is tracked independently, no module-level mutable state, all working
  storage local or in per-path objects. (Stateful user evaluator types
  will then need one copy per thread.)
- **No odeint dependency** (checked 2026-09-30): its RHS interface
  `dydt(y, t)` is real-valued and has no way to pass context (the
  homotopy object) except module-level state, and path tracking needs
  a corrector-driven step control, not ODE error control. The
  predictor is implemented directly in the tracker. Reconsider if
  odeint gets a complex / object-based (class with deferred `dydt`)
  interface.
- **Metadata**: MIT license, author "Thomas Biekötter", maintainer
  thomas.biekoetter@desy.de (as in fortransitions).
- **README**: states that the code was written by Claude (Anthropic),
  like `fortransitions/README.md`.

## Status and next steps

**Milestone 1 (done, 2026-09-30)**: total-degree homotopy solver for
regular solutions: all modules of the architecture table, example in
`app/example.f90`, tests for every module plus Katsura-3/4 against
sympy references. Known limitation: without an end game, paths to
singular endpoints (multiple roots, singular points at infinity) fail
close to t = 0 (`path_failed_min_step`, t ~ 1e-14 .. 1e-7) or end
flagged singular; they are excluded from `solutions`.

**Milestone 2 (next)**: end game near t = 0 (power-series or Cauchy end
game, winding numbers), so that singular endpoints and points at
infinity are classified properly; clustering of singular solutions with
multiplicities (`singular_solutions`); then consider re-tracking paths
with tighter settings when `nduplicates > 0`. Later items: see the
roadmap above (scaling, polyhedral homotopy, parameter homotopy,
monodromy).

## Working conventions

- Keep the README up to date with what is implemented and how to use it
  (see `fortransitions/README.md` for the expected level of detail).
- **Commit in fine-grained steps**: one commit per self-contained unit
  of work (e.g. project setup; one new module together with its test;
  one feature or fix), each leaving the project building and all tests
  passing with gfortran (debug and release profiles; ifx too when
  available). Never commit a broken build. Descriptive commit messages;
  git identity is already configured. Do not push unless asked.
