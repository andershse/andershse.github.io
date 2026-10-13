Generic Local Identifiability Does Not Imply Generic Global Identifiability
in Linear Networks
Draft for IEEE Control Systems Letters
Revised following the supplied referee report, October 8, 2026

Main files
----------
local_global_identifiability.tex         Manuscript, including both diagrams
local_global_identifiability_submission.pdf Current manuscript without a pending archive citation
ieeeconf.cls                            Unmodified official PaperPlaza class
REVISION_NOTES.txt                      Response to the referee's suggestions
identifiability-verification.zip         Prepared public verification archive
ZENODO_SETUP.txt                        Account and publication instructions
local_global_identifiability_source.zip Manuscript source package

Earlier PDFs have been retained. The current manuscript is the _submission PDF.
The pre-revision manuscript and README are in pre_review_20261008/.

Build
-----
Run twice from the directory containing the manuscript:

  pdflatex -jobname=local_global_identifiability_submission -interaction=nonstopmode -halt-on-error local_global_identifiability.tex

The bibliography is included in the LaTeX file; BibTeX is not required.
Standard packages include TikZ, balance, amsmath, amssymb, and hyperref.
The document uses the official 10-point US Letter class and margin commands.
No custom margins, font sizes, or line spacing were used to meet the limit.
If verification_doi.tex is present, its DOI is included automatically.

Template and journal instructions, checked October 8, 2026:
https://css.paperplaza.net/conferences/support/files/ieeeconf.zip
https://ieeecss.org/publication/ieee-control-systems-letters/author-information

Author information
------------------
Anders Hansson, Senior Member, IEEE
Division of Automatic Control, Linköping University
Email: anders.g.hansson@liu.se

The manuscript has not been submitted. The acknowledgment describes the
AI assistance; a short citation in each section's footnote implements the
journal's section-level disclosure requirement. Section headings are clean.

Verification
------------
The prepared archive contains both verification scripts, pinned dependencies,
exact JSON outputs, README, citation metadata, MIT license, and a SHA-256
manifest. It excludes the manuscript, referee report, and account data.
From the extracted archive directory:

  python -m pip install -r requirements.txt
  python verify_counterexample.py
  python verify_spectral_counterexample.py

In the workspace, the current scripts and outputs are in ../research/.
Identical copies are in verification/ and the verification ZIP.
Both scripts passed during this revision. Checks include the symbolic
involution and invariant domain, nonsymmetric parameters, equality of the
transfer matrices, nilpotence, exact power bounds, and full Jacobian column
rank 42 at both witnesses. Derivatives use all independent edge positions.
The spectral script also verifies a perturbation with 42 distinct edge
values, the exact example with only one admissible solution, and the
conditioning formula. See verification/README.txt for the certificates.
All mathematical proofs needed for the paper are in the manuscript.

Public archive status
---------------------
The Zenodo payload, metadata, and uploader are prepared. No public record
or DOI has yet been created. Publication requires authentication through
the author's Zenodo account; see ZENODO_SETUP.txt. The uploader was tested
with simulated API responses and an offline archive check. A live upload
has not been performed. A public code archive and DOI are optional for
initial L-CSS submission. The manuscript currently states that code,
dependencies, and exact outputs are available from the author on request.
It includes the software citation only when a public DOI is supplied.
The self-contained manuscript can be submitted without a Zenodo deposit.

Scope and editorial audit
-------------------------
The paper was audited against EDITORIAL_RULES.md: definitions and dependency
order, consistency of symbols, terminology, and the scope of each claim.
The scalar result uses independent complex coefficients with volume in
2m real coordinates; Theorem 1 also states the real-domain conclusion.
The theorem disproves generic global identifiability under a common
spectral-radius restriction. It does not claim that almost every admissible
network is ambiguous: Section V-B proves open regions with one solution.

Section VI treats whole rational functions as coordinates over C(lambda).
The generic algebraic fiber has exactly two points, and the rational
transformation preserves rational-function-valued coordinates. A real open
set of constant witnesses is Zariski dense over the algebraic closure of
C(lambda), as is its image in the transfer-matrix image. This rules out
algebraically generic uniqueness under properness, stability, and the
spectral-radius bound. Genericity in prescribed positive-order numerator
and denominator coefficient families and strictly proper networks remain
separate questions. The construction supplies no uniform conditioning
bound or neighborhood size as the spectral-radius bound decreases.

The manuscript uses "transfer matrix" for C(I-G)^(-1)B, F, and H and
relates these to the cited algebraic paper's closed-loop transfer function.
The counts of six and twelve nodes are upper bounds, not minimality claims.

Format verification
-------------------
The abstract has 174 whitespace-separated words. The revised PDF has six
US Letter pages, including two diagrams and references. All fonts are
embedded Type 1 fonts. All six pages were visually checked.
The final LaTeX log has no undefined references, overfull or underfull boxes,
or package warnings. MiKTeX prints an environment maintenance notice that
an administrator has not checked for updates; compilation succeeds.
The final columns are balanced. Adding the public DOI requires rebuilding
and checking the PDF again.
