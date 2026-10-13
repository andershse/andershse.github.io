r"""Exact counterexample to local => global with rho(G) < 1.

The model has 12 nodes, zero diagonal, 42 independent unknown nonzero
entries, 9 excited nodes, and 9 measured nodes. This verifies a rational
pair, strict spectral-radius bounds, full local rank, and the augmentation
identity. Run with the same environment as verify_counterexample.py.
"""

import json
from pathlib import Path

import numpy as np
import sympy as s


def mod_matrix(M, prime):
    return [[int(s.numer(t)) % prime * pow(int(s.denom(t)) % prime, -1, prime) % prime
             for t in row] for row in M.tolist()]


def rank_mod(rows, prime):
    rows = [row[:] for row in rows]
    n, m = len(rows), len(rows[0])
    rank = 0
    for col in range(m):
        pivot = next((i for i in range(rank,n) if rows[i][col] % prime), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        scale = pow(rows[rank][col], -1, prime)
        rows[rank] = [(v*scale) % prime for v in rows[rank]]
        for i in range(rank+1,n):
            coefficient = rows[i][col]
            if coefficient:
                rows[i] = [(v-coefficient*w) % prime
                           for v,w in zip(rows[i],rows[rank])]
        rank += 1
        if rank == n:
            break
    return rank


U = s.Matrix([[0,0,1],[1,0,0],[0,1,0]])
I = s.eye(3)
Z = s.zeros(3)
S = U + U.T
B0 = I.col_join(Z)
C0 = Z.row_join(I)
Kc = S.row_join(-U)
Lc = U.T.col_join(-S)
Jc = U.T.row_join(-S).col_join(I.row_join(U))


def build(r, x):
    original = (-r*U).row_join(Z).col_join((x*I).row_join(-r*U.T))
    augmented = original.row_join(B0*Kc).col_join((Lc*C0).row_join(Jc))
    return original, augmented


B = s.eye(12)[:, [0,1,2,6,7,8,9,10,11]]
C = s.eye(12)[[3,4,5,6,7,8,9,10,11], :]
r = s.Rational(99,100)
original1, G1 = build(r, s.Integer(1))
original2, G2 = build(1/r, 1/r**2)
_, G0 = build(s.Integer(1), s.Integer(1))
assert G0**4 == s.zeros(12)
assert G0**3 != s.zeros(12)
assert max(sum(abs(G0[i,j]) for j in range(12)) for i in range(12)) == 4
assert G1 != G2
edges = [(i,j) for i in range(12) for j in range(12) if G1[i,j] != 0]
assert len(edges) == 42
assert all(G1[i,i] == G2[i,i] == 0 for i in range(12))
assert edges == [(i,j) for i in range(12) for j in range(12) if G2[i,j] != 0]
print("Same zero-diagonal pattern, 42 independent edge positions: PASS",flush=True)
print("Coalescence matrix G0^4 = 0: PASS",flush=True)

F = C0 * (s.eye(6)-original1).inv() * B0
assert F == C0 * (s.eye(6)-original2).inv() * B0
records = []
observations = []
prime = 1000003
for label,G in [("first",G1),("second",G2)]:
    power = G**4
    norm4 = max(sum(abs(power[i,j]) for j in range(12)) for i in range(12))
    assert norm4 < 1
    T = (s.eye(12)-G).inv()
    H = C*T*B
    observations.append(H)
    inverse_expected = F.inv().row_join(-Kc).col_join((-Lc).row_join(s.eye(6)-Jc))
    assert H.inv() == inverse_expected
    # Exact modular certificate: a nonzero minor modulo this prime is
    # a nonzero rational minor, proving full column rank over Q and C.
    CT = mod_matrix(C*T, prime)
    TB = mod_matrix(T*B, prime)
    jacobian = [[CT[c][i]*TB[j][b] % prime for i,j in edges]
                for c in range(9) for b in range(9)]
    jacobian_rank = rank_mod(jacobian, prime)
    assert jacobian_rank == 42
    approx_rho = float(max(abs(np.linalg.eigvals(np.array(G).astype(float)))))
    record = dict(label=label, norm_fourth_power=str(norm4),
                  norm_fourth_power_decimal=float(norm4),
                  spectral_radius_numeric=approx_rho,
                  certified_spectral_radius_upper_bound=float(norm4)**0.25,
                  jacobian_rank_mod_prime=jacobian_rank,
                  jacobian_prime=prime,
                  determinant_I_minus_G=str((s.eye(12)-G).det()),
                  matrix=[[str(e) for e in row] for row in G.tolist()])
    records.append(record)
    print(json.dumps({k:v for k,v in record.items() if k!='matrix'},indent=2),flush=True)

assert observations[0] == observations[1]
print("Exact equality of full 9 by 9 observations: PASS",flush=True)

# A distinct, freely perturbed point checks the same mechanism without
# any of the numerical coincidences in the compact examples.
aa = [-r, -r-s.Rational(1,100000), -r+s.Rational(2,100000)]
dd = [-r+s.Rational(3,100000), -r-s.Rational(4,100000), -r+s.Rational(5,100000)]
xx = [1+s.Rational(6,100000), 1-s.Rational(7,100000), 1+s.Rational(8,100000)]
k = [aa[i]*dd[i]*xx[(i+1)%3]/xx[i] for i in range(3)]
h = [1+k[i]+k[i]*k[(i+1)%3] for i in range(3)]
ratios = [s.cancel(k[i]*h[(i+1)%3]/h[i]) for i in range(3)]

def full(a,d,x,ctrlK,ctrlL,ctrlJ):
    A, D = s.zeros(3),s.zeros(3)
    for i in range(3):
        A[(i+1)%3,i] = a[i]
        D[i,(i+1)%3] = d[i]
    O=A.row_join(Z).col_join(s.diag(*x).row_join(D))
    return O.row_join(B0*ctrlK).col_join((ctrlL*C0).row_join(ctrlJ))

controller_matrices = [Kc.copy(),Lc.copy(),Jc.copy()]
counter = 10
for matrix in controller_matrices:
    for i in range(matrix.rows):
        for j in range(matrix.cols):
            if matrix[i,j] != 0:
                matrix[i,j] += s.Rational(counter,100000)
                counter += 1
perturbed1 = full(aa,dd,xx,*controller_matrices)
perturbed2 = full([aa[i]/ratios[i] for i in range(3)],
                 [dd[i]/ratios[i] for i in range(3)],
                 [xx[i]/ratios[i] for i in range(3)],*controller_matrices)
for M in [perturbed1,perturbed2]:
    MP = M**4
    assert max(sum(abs(MP[i,j]) for j in range(12)) for i in range(12)) < 1
assert C*(s.eye(12)-perturbed1).inv()*B == C*(s.eye(12)-perturbed2).inv()*B
assert len(set(perturbed1[i,j] for i,j in edges)) == 42
print("Rational perturbation with 42 distinct nonzero edge values: PASS",flush=True)

# A full-pattern witness with exactly one spectrally admissible candidate.
# The base two-solution theorem and augmentation lemma exclude other candidates.
small = s.Rational(1, 1000)
unique1 = full([-s.Rational(1,2)]*3, [-s.Rational(1,2)]*3,
               [s.Integer(1)]*3, small*Kc, small*Lc, small*Jc)
unique2 = full([-s.Integer(2)]*3, [-s.Integer(2)]*3,
               [s.Integer(4)]*3, small*Kc, small*Lc, small*Jc)
for M in [unique1, unique2]:
    assert [(i,j) for i in range(12) for j in range(12) if M[i,j] != 0] == edges
assert C*(s.eye(12)-unique1).inv()*B == C*(s.eye(12)-unique2).inv()*B
unique_power = unique1**4
unique_norm = max(sum(abs(unique_power[i,j]) for j in range(12))
                  for i in range(12))
excluded_trace_bound = abs(s.trace(unique2**3))/12
assert unique_norm == s.Rational(564741, 1000000) < 1
assert excluded_trace_bound == s.Rational(3999994001, 1000000000) > 1
unique_record = dict(admissible_fourth_power_norm=str(unique_norm),
                     excluded_spectral_radius_cubed_lower_bound=str(excluded_trace_bound),
                     matrices=[[[str(e) for e in row] for row in M.tolist()]
                               for M in [unique1, unique2]])
print("One admissible candidate: exact power and trace certificates PASS", flush=True)

# Along the symmetric family, the reduced derivative has one vanishing
# singular value as t approaches 1. Compute its Gram-matrix spectrum exactly.
ts = s.symbols("t", positive=True)
reduced = I - ts**2*U
gram_eigenvalues = (reduced.T*reduced).eigenvals()
assert gram_eigenvalues == {
    ts**4 - 2*ts**2 + 1: 1,
    ts**4 + ts**2 + 1: 2,
}
print("Reduced derivative singular values and limiting infinity norm: PASS", flush=True)

result = dict(nodes=12,independent_edges=42,excited_nodes=[1,2,3,7,8,9,10,11,12],
              measured_nodes=[4,5,6,7,8,9,10,11,12],same_observation_exact=True,
              inverse_observation=[[str(e) for e in row] for row in observations[0].inv().tolist()],
              examples=records, distinct_edge_perturbation_verified=True,
              unique_admissible_example=unique_record,
              limiting_matrix_infinity_norm=4, conditioning_identity_verified=True)
Path(__file__).with_name("spectral_verification_results.json").write_text(json.dumps(result,indent=2))
print("ALL EXACT CERTIFICATES PASSED",flush=True)
