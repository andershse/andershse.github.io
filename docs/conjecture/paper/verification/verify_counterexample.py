r"""Exact certificates for a six-node local/global identifiability counterexample.

Run: .venv\Scripts\python.exe research\verify_counterexample.py
All nine edges in the model are independent unknowns. Symmetric numerical
specializations below are examples only, not constraints on the model.
"""

import json
from pathlib import Path

import sympy as s


def blocks(a, d, x):
    A = s.zeros(3)
    D = s.zeros(3)
    for i in range(3):
        A[(i + 1) % 3, i] = a[i]
        D[i, (i + 1) % 3] = d[i]
    X = s.diag(*x)
    G = A.row_join(s.zeros(3)).col_join(X.row_join(D))
    B = s.eye(3).col_join(s.zeros(3))
    C = s.zeros(3).row_join(s.eye(3))
    return A, D, X, G, B, C


def zero_matrix(M):
    return all(s.cancel(e) == 0 for e in M)


a = s.symbols("a1:4", nonzero=True)
d = s.symbols("d1:4", nonzero=True)
x = s.symbols("x1:4", nonzero=True)
A, D, X, G, B, C = blocks(a, d, x)
I = s.eye(3)
K = (I - A) * X.inv() * (I - D)
print("Inverse observed transfer matrix K:", K, flush=True)

# Coordinate change u_i=a_i/x_i, v_i=d_i/x_i, z_i=1/x_i.
u = s.symbols("u1:4")
v = s.symbols("v1:4")
z = s.symbols("z1:4", nonzero=True)
diag = s.Matrix([z[i] + u[(i - 1) % 3] * v[(i - 1) % 3] / z[(i - 1) % 3]
                 for i in range(3)])
det_reduced = s.factor(diag.jacobian(z).det())
assert s.cancel(det_reduced - (1 - s.prod(u) * s.prod(v) / s.prod(z)**2)) == 0
print("Reduced Jacobian determinant:", det_reduced, flush=True)

# Rational formula for the second point of the fiber.
k = [a[i] * d[i] * x[(i + 1) % 3] / x[i] for i in range(3)]
h = [1 + k[i] + k[i] * k[(i + 1) % 3] for i in range(3)]
r = [s.cancel(k[i] * h[(i + 1) % 3] / h[i]) for i in range(3)]
ap = [a[i] / r[i] for i in range(3)]
dp = [d[i] / r[i] for i in range(3)]
xp = [x[i] / r[i] for i in range(3)]
Ap, Dp, Xp, Gp, _, _ = blocks(ap, dp, xp)
Kp = (I - Ap) * Xp.inv() * (I - Dp)
assert zero_matrix(Kp - K)
assert s.cancel(s.prod(ap) - 1 / s.prod(d)) == 0
assert s.cancel(s.prod(dp) - 1 / s.prod(a)) == 0
print("Symbolic equality of inverse transfer matrices: PASS", flush=True)
print("Cycle products exchange with reciprocals: PASS", flush=True)

# Check invariance of the stated domain and the involution in auxiliary
# coordinates. Nonzero k_i and h_i imply nonzero transformed h_i.
ks = s.symbols("k1:4", nonzero=True)
hs = [1 + ks[i] + ks[i] * ks[(i + 1) % 3] for i in range(3)]
etas = [ks[i] * hs[(i + 1) % 3] / hs[i] for i in range(3)]
kps = [s.cancel(ks[i] / (etas[i] * etas[(i + 1) % 3]))
       for i in range(3)]
hps = [s.cancel(1 + kps[i] + kps[i] * kps[(i + 1) % 3])
       for i in range(3)]
for i in range(3):
    assert s.cancel(hps[i] - hs[(i + 1) % 3] /
                    (ks[(i + 1) % 3] * ks[(i + 2) % 3])) == 0
    etap = s.cancel(kps[i] * hps[(i + 1) % 3] / hps[i])
    assert s.cancel(etap * etas[i] - 1) == 0
    assert s.cancel(r[i] - 1 - (s.prod(a) * s.prod(d) - 1) / h[i]) == 0
kc = [u[i] * v[i] / (z[i] * z[(i + 1) % 3]) for i in range(3)]
leading = u[1] * v[1] - diag[1] * diag[2]
assert s.cancel(leading + z[1] * z[2] * (1 + kc[0] + kc[0] * kc[1])) == 0
print("Invariant domain, rational involution, and quadratic leading coefficient: PASS",
      flush=True)

# The fixed-point equation for the three successive Mobius transformations
# has degree two. Verify with independent observed entries.
e1, e2, e3, b1, b2, b3, t = s.symbols("e1 e2 e3 b1 b2 b3 t")
z2 = e2 - b1 / t
z3 = e3 - b2 / z2
fixed = s.factor(s.together(t - (e1 - b3 / z3)))
num, den = s.fraction(fixed)
assert s.degree(num, t) == 2
print("Generic quadratic:", s.Poly(num, t).as_expr(), flush=True)

# Independent exact checks at a simple point and a point with distinct edges.
records = []
for aa, dd, xx in [
    ([s.Rational(1, 2)] * 3, [s.Rational(1, 3)] * 3, [s.Integer(1)] * 3),
    ([s.Rational(1, 2), s.Rational(2, 3), s.Rational(3, 4)],
     [s.Rational(1, 5), s.Rational(2, 7), s.Rational(3, 11)],
     [s.Rational(2), s.Rational(3), s.Rational(5)]),
]:
    sub = dict(zip(a + d + x, aa + dd + xx))
    GG = G.subs(sub)
    GGp = Gp.subs(sub).applyfunc(s.cancel)
    assert GG != GGp
    assert all((GG[i, j] == 0) == (GGp[i, j] == 0)
               for i in range(6) for j in range(6))
    T = (s.eye(6) - GG).inv()
    Tp = (s.eye(6) - GGp).inv()
    F = C * T * B
    assert C * Tp * B == F
    edges = [(i, j) for i in range(6) for j in range(6) if G[i, j] != 0]
    J = s.Matrix.hstack(*[
        s.Matrix(C * T[:, i] * T[j, :] * B).reshape(9, 1)
        for i, j in edges
    ])
    Jp = s.Matrix.hstack(*[
        s.Matrix(C * Tp[:, i] * Tp[j, :] * B).reshape(9, 1)
        for i, j in edges
    ])
    assert J.det() != 0 and Jp.det() != 0
    rec = {
        "original_parameters": [str(q) for q in aa + dd + xx],
        "alternative_parameters": [str(s.cancel(q.subs(sub))) for q in ap + dp + xp],
        "G": str(GG),
        "G_alternative": str(GGp),
        "observed_transfer": str(F),
        "det_I_minus_G": str((s.eye(6) - GG).det()),
        "det_I_minus_G_alternative": str((s.eye(6) - GGp).det()),
        "jacobian_det": str(J.det()),
        "alternative_jacobian_det": str(Jp.det()),
    }
    records.append(rec)
    print(json.dumps(rec, indent=2), flush=True)

# A nonconstant specialization with proper, stable edges and internal resolvents.
# Its poles lie among 0, -1/4, -1/5, and roots of the two cycle polynomials.
# The note bounds those cycle roots by 1/2 and 3/10, respectively.
q = s.symbols("q")
aq = q / (2 * (q + s.Rational(1, 4)))
dq = q / (3 * (q + s.Rational(1, 5)))
Aq, Dq, Xq, _, _, _ = blocks([aq] * 3, [dq] * 3, [s.Integer(1)] * 3)
Aqp, Dqp, Xqp, _, _, _ = blocks([1/dq] * 3, [1/aq] * 3, [1/(aq*dq)] * 3)
TA, TD = (I-Aq).inv(), (I-Dq).inv()
TAp, TDp = (I-Aqp).inv(), (I-Dqp).inv()
Fq, Fqp = TD * Xq * TA, TDp * Xqp * TAp
assert zero_matrix(Fq - Fqp)
stable_polynomial = s.Poly(
    q**3 * (q+s.Rational(1,4))**3 * (q+s.Rational(1,5))**3
    * (8*(q+s.Rational(1,4))**3-q**3)**3
    * (27*(q+s.Rational(1,5))**3-q**3)**3, q)
for M in [Aq, Dq, Xq, Aqp, Dqp, Xqp, TA, TD, TAp, TDp, Fq]:
    for entry in M:
        numerator, denominator = s.fraction(s.cancel(entry))
        assert s.degree(numerator, q) <= s.degree(denominator, q)
        assert stable_polynomial.rem(s.Poly(denominator, q)).is_zero
print("Nonconstant rational example: equality, properness and stable pole set PASS", flush=True)

Path(__file__).with_name("verification_results.json").write_text(
    json.dumps({"symbolic_identity": True, "fiber_degree": 2,
                "invariant_domain_verified": True, "involution_verified": True,
                "leading_coefficient_identity_verified": True,
                "examples": records, "nonconstant_stable_example_verified": True}, indent=2)
)
print("ALL EXACT CHECKS PASSED", flush=True)
