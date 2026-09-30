# Reference scripts

These Python scripts compute reference solutions for the Fortran test
programs independently of homotopy continuation (Groebner bases with
sympy, high-precision root finding with mpmath), and print them as
Fortran array constructors that are pasted into the tests:

| Script           | Produces reference values for |
| ---------------- | ----------------------------- |
| `ref_katsura.py` | `test/katsura.f90`            |

Run with

```sh
python3 ref_katsura.py 3 4
```

The values in the tests were generated with Python 3.13.11,
sympy 1.14.0 and mpmath 1.3.0 (installed in the pyenv virtualenv
`polyhomocont`). The solutions are computed with 40 significant digits
and printed with 17, so they are exact to double precision; the tests
compare solution sets up to ordering (see `set_distance` in
`test/testing.f90`).
