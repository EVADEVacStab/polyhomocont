"""Reference solutions of the Katsura-n systems for test/katsura.f90.

Katsura-n has the n + 1 unknowns u_0, ..., u_n and the equations

    u_0 + 2 (u_1 + ... + u_n) - 1 = 0,
    sum_{l=-n}^{n} u_{|l|} u_{|m-l|} - u_m = 0,   m = 0, ..., n - 1,

with u_k = 0 for k > n. It has 2^n solutions, all finite and regular,
which equals its Bezout number (1 * 2^n).

The solutions are computed independently of homotopy continuation: a
lexicographic Groebner basis (sympy) in shape position,
    u_i - g_i(u_n) (i < n),   p(u_n),
the roots of the univariate polynomial p with mpmath at 40 significant
digits, and back-substitution. The script prints the solutions as
Fortran array constructors, pasted into test/katsura.f90.
"""

import sys

import mpmath
import sympy

mpmath.mp.dps = 40


def katsura(n):
    u = sympy.symbols(f"u0:{n + 1}")

    def uu(k):
        k = abs(k)
        return u[k] if k <= n else 0

    eqs = [u[0] + 2 * sum(u[1:]) - 1]
    for m in range(n):
        eqs.append(sympy.expand(
            sum(uu(l) * uu(m - l) for l in range(-n, n + 1)) - u[m]))
    return u, eqs


def solve_shape_position(u, eqs):
    """Solutions of a zero-dimensional system whose lex Groebner basis
    (u_0 > ... > u_n) is in shape position."""
    gb = sympy.groebner(eqs, *u, order="lex")
    polys = list(gb.exprs)
    last = u[-1]
    univariate = [p for p in polys if p.free_symbols <= {last}]
    assert len(univariate) == 1, "basis not in shape position"
    p = sympy.Poly(univariate[0], last)
    subs = {}
    for var in u[:-1]:
        cand = [q for q in polys if sympy.Poly(q, var).degree() == 1
                and q.free_symbols <= {var, last}]
        assert cand, "basis not in shape position"
        g = sympy.solve(cand[0], var)[0]
        subs[var] = sympy.lambdify(last, g, "mpmath")
    coeffs = [mpmath.mpmathify(sympy.N(c, 60)) for c in p.all_coeffs()]
    roots = mpmath.polyroots(coeffs, maxsteps=500, extraprec=200)
    sols = []
    for r in roots:
        sols.append([subs[var](r) for var in u[:-1]] + [r])
    return sols, p.degree()


def residual(u, eqs, sol):
    f = [sympy.lambdify(u, e, "mpmath") for e in eqs]
    return max(abs(fi(*sol)) for fi in f)


def fortran_complex(z):
    return f"({float(z.real):.17e}_wp, {float(z.imag):.17e}_wp)"


def main():
    ns = [int(a) for a in sys.argv[1:]] or [3, 4]
    for n in ns:
        u, eqs = katsura(n)
        sols, deg = solve_shape_position(u, eqs)
        assert deg == 2**n and len(sols) == 2**n
        res = max(residual(u, eqs, s) for s in sols)
        nreal = sum(all(abs(mpmath.mpc(c).imag) < mpmath.mpf("1e-30")
                        for c in s) for s in sols)
        print(f"! Katsura-{n}: {len(sols)} solutions ({nreal} real), "
              f"max residual {float(res):.1e}")
        print(f"sols{n} = reshape([  &")
        entries = [fortran_complex(mpmath.mpc(c)) for s in sols for c in s]
        print(",  &\n".join(f"  {e}" for e in entries)
              + f"], [{n + 1}, {2**n}])")


if __name__ == "__main__":
    main()
