"""Search DAGs for a second isolated solution of the same observation map.

Independent unknown edges only, zero diagonal, fixed node selections.
The block triangular network has nilpotent A and D, so rho(G)=0 exactly.
Numerical findings must subsequently be certified symbolically.
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.optimize import root


def run(p, trials, starts, seed):
    rng = np.random.default_rng(seed)
    possible = [(i,j) for i in range(p) for j in range(i)]
    count_local = 0
    beginning = time.monotonic()
    for trial in range(trials):
        na = int(rng.integers(2, len(possible)+1))
        nd = int(rng.integers(2, len(possible)+1))
        ae = [possible[i] for i in rng.choice(len(possible), na, replace=False)]
        de = [possible[i] for i in rng.choice(len(possible), nd, replace=False)]
        zeros = [divmod(int(k),p) for k in rng.choice(p*p, na+nd, replace=False)]
        zr, zc = zip(*zeros)
        zr, zc = np.array(zr), np.array(zc)
        theta = rng.integers(1,8,na+nd).astype(float)
        A, D = np.zeros((p,p)), np.zeros((p,p))
        for value, (i,j) in zip(theta[:na], ae): A[i,j] = value
        for value, (i,j) in zip(theta[na:], de): D[i,j] = value
        X = rng.integers(1,8,(p,p)).astype(float)
        X[zr, zc] = 0
        F = np.linalg.solve(np.eye(p)-D, X) @ np.linalg.inv(np.eye(p)-A)

        def fun(t):
            AA, DD = np.zeros((p,p)), np.zeros((p,p))
            for value, (i,j) in zip(t[:na], ae): AA[i,j] = value
            for value, (i,j) in zip(t[na:], de): DD[i,j] = value
            return ((np.eye(p)-DD) @ F @ (np.eye(p)-AA))[zr,zc]

        def jac(t):
            AA, DD = np.zeros((p,p)), np.zeros((p,p))
            for value, (i,j) in zip(t[:na], ae): AA[i,j] = value
            for value, (i,j) in zip(t[na:], de): DD[i,j] = value
            LF = (np.eye(p)-DD) @ F
            FR = F @ (np.eye(p)-AA)
            result = np.empty((na+nd,na+nd))
            for k,(i,j) in enumerate(ae): result[:,k] = -LF[zr,i]*(zc==j)
            for k,(i,j) in enumerate(de): result[:,na+k] = (zr==i)*(-FR[j,zc])
            return result

        sv = np.linalg.svd(jac(theta), compute_uv=False)
        if sv[-1] < 1e-9*sv[0]: continue
        count_local += 1
        for start in range(starts):
            initial = rng.normal(0, 5 if start%2 else 15, na+nd)
            sol = root(fun, initial, jac=jac, options={'xtol':1e-10})
            if (np.linalg.norm(fun(sol.x)) < 1e-7 and np.linalg.norm(sol.x-theta)>1e-4
                and np.max(np.abs(sol.x))<1e5
                and np.linalg.matrix_rank(jac(sol.x),tol=1e-7)==na+nd):
                result = dict(p=p, trial=trial, local_count=count_local,
                              a_edges=ae,d_edges=de,x_zeros=zeros,
                              theta=theta.tolist(),alternative=sol.x.tolist(),
                              X=X.tolist(),F=F.tolist(),residual=float(np.linalg.norm(fun(sol.x))))
                path=Path(__file__).with_name(f'dag_candidate_{p}.json')
                path.write_text(json.dumps(result,indent=2))
                print('CANDIDATE '+json.dumps(result),flush=True)
                return
        if count_local%25==0:
            print(f'p={p} trials={trial+1} local={count_local} seconds={time.monotonic()-beginning:.1f}',flush=True)
    print(f'FINISHED p={p} trials={trials} local={count_local} seconds={time.monotonic()-beginning:.1f}',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--p',type=int,default=4)
    parser.add_argument('--trials',type=int,default=10000)
    parser.add_argument('--starts',type=int,default=15)
    parser.add_argument('--seed',type=int,default=20261008)
    args=parser.parse_args()
    run(args.p,args.trials,args.starts,args.seed)
