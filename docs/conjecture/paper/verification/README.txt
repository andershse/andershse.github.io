Verification code for generic local and global identifiability in linear networks
Version 1.0.0
Author: Anders Hansson
Division of Automatic Control, Linköping University
Contact: anders.g.hansson@liu.se

Companion manuscript
--------------------
Generic Local Identifiability Does Not Imply Generic Global Identifiability
in Linear Networks

This archive contains the arithmetic checks for the manuscript. Its main
proofs are in the manuscript. The checks distinguish the full independent
edge parameterization from the symmetric numerical witness family.

Run
---
Use Python 3.11 or later. From the extracted archive directory:

  python -m pip install -r requirements.txt
  python verify_counterexample.py
  python verify_spectral_counterexample.py

Both scripts must end with a message that all exact checks passed.
They rewrite verification_results.json and spectral_verification_results.json
in their own directory. The supplied JSON files are outputs from a successful
run with the pinned dependencies.

Checks
------
verify_counterexample.py:
  - Reconstructs the six-node graph and selection matrices.
  - Verifies the equality of the two inverse transfer matrices symbolically.
  - Verifies the rank-nine determinant and reciprocal cycle products.
  - Checks the invariant domain, rational involution, and leading coefficient
    of the quadratic giving the fiber bound.
  - Checks symmetric and nonsymmetric rational-number examples.
  - Includes an additional nonconstant proper example with stable pole sets.
    This last example does not assert that both six-node network matrices
    satisfy the spectral-radius bound.

verify_spectral_counterexample.py:
  - Reconstructs all 42 permitted edge positions and the selections.
  - Verifies that the limiting matrix has fourth power zero, third power
    nonzero, and infinity norm four.
  - Checks exact equality of the full nine-by-nine transfer matrices.
  - Computes exact fourth-power norm certificates.
  - Certifies Jacobian column rank 42 at both rational witnesses.
    For each permitted entry G[i,j], it forms the full derivative
    C(I-G)^(-1)[:,i] (I-G)^(-1)[j,:]B. These are derivatives with
    respect to all independent edges, not derivatives along the parameter t.
    Rank is computed exactly modulo the prime 1000003. A nonzero minor
    modulo this prime certifies a nonzero rational minor.
  - Checks a rational perturbation with 42 distinct edge values.
  - Checks a second exact pair having just one spectrally admissible member.
  - Verifies the Gram-matrix spectrum giving the conditioning formula.

Expected certificates
---------------------
Two admissible networks:
  ||G1^4||_infinity = 33910001/100000000
  ||G2^4||_infinity = 16976126/32019867
  Both full Jacobian column ranks are 42.

One admissible network and one excluded alternative:
  ||G1^4||_infinity = 564741/1000000 < 1
  |trace(G2^3)|/12  = 3999994001/1000000000 > 1
  The second quantity is a lower bound for rho(G2)^3.

Numerical eigenvalues printed by NumPy are supplementary. The admissibility
and rank certificates use exact arithmetic.

Scope
-----
The scalar theorem concerns independent real or complex edge coefficients.
The field extension treats whole rational transfer functions as coordinates.
The archive does not establish genericity in every prescribed numerator and
denominator coefficient family or the case of strictly proper modules.

Files and integrity
-------------------
The archive contains this README, two Python scripts, requirements.txt,
their two JSON outputs, CITATION.cff, LICENSE, and MANIFEST.json.
MANIFEST.json records the SHA-256 hash of each other file in the archive.
No manuscript source, referee correspondence, credentials, or browser data
are included.

License and assistance
----------------------
The verification code is prepared for release under the MIT license.
The dependencies have their own licenses.
OpenAI Codex assisted in developing the construction, proofs, verification
scripts, and documentation.
