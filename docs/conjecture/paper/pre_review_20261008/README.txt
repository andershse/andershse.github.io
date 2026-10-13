Local and Global Identifiability Can Differ in Linear Networks
Draft for IEEE Control Systems Letters

Main files
----------
local_global_identifiability.tex : manuscript and embedded vector diagram
local_global_identifiability_algebraic.pdf : current compiled manuscript
ieeeconf.cls                    : unmodified official PaperPlaza class

The current PDF uses the _algebraic filename because Windows had the
original local_global_identifiability.pdf open during the update.
The LaTeX source retains its original filename.

Build
-----
Run these commands from the directory containing the manuscript:

  pdflatex -jobname=local_global_identifiability_algebraic -interaction=nonstopmode -halt-on-error local_global_identifiability.tex
  pdflatex -jobname=local_global_identifiability_algebraic -interaction=nonstopmode -halt-on-error local_global_identifiability.tex

The bibliography is included in the LaTeX file; BibTeX is not required.
The document uses standard LaTeX packages, including TikZ and balance.
The class, font size, paper size, and margin commands follow the official
template. No custom changes to margins, font sizes, or line spacing were
used to meet the page limit.

Template and journal requirements, checked October 8, 2026:
https://css.paperplaza.net/conferences/support/files/ieeeconf.zip
https://ieeecss.org/publication/ieee-control-systems-letters/author-information

Author information
------------------
Anders Hansson, Senior Member, IEEE
Division of Automatic Control, Linköping University
Email: anders.g.hansson@liu.se

The manuscript has not been submitted.

The acknowledgment describes the AI assistance. The citations in section
headings implement the journal's requirement for a citation in each section
using AI-generated text.

Verification
------------
In the workspace, the following files are in ../research/. They are also
included at the top level of the source archive:

  verify_counterexample.py
  verify_spectral_counterexample.py
  requirements.txt

Install the dependencies with:

  python -m pip install -r requirements.txt

Then run:

  python verify_counterexample.py
  python verify_spectral_counterexample.py

The scripts check the symbolic alternative, the reduced Jacobian
determinant, the quadratic equation, exact equality of transfer matrices,
the matrix-power bounds, and the full Jacobian ranks at the rational pair.
The second script also checks a perturbation with 42 distinct edge values.
Each script writes its verification results as a JSON file.
Both scripts passed during preparation of this draft.

Editorial and format review
---------------------------
The manuscript was reviewed against EDITORIAL_RULES.md, including the
order of definitions, the roles of symbols, consistency of notation, and
the distinction between scalar and algebraic genericity.
The scalar parameter domain is C^m, matching the cited
fixed-frequency formulation. Volume is measured in the 2m real and
imaginary coordinates, and the theorem also holds on the real parameter
domain. The explicit matrices have real rational entries. The added node
signals are denoted w_c to avoid collision with the auxiliary scalars h_i.

Section VI and Proposition 2 extend the two-solution result to the field
C(lambda), with lambda a symbolic frequency variable. The generic
algebraic fiber has exactly two points, and both points have rational
transfer-function coordinates whenever one does. The proof also excludes
additional augmented solutions outside the block-inverse domain.
The real open set of constant examples is Zariski dense over the
algebraic closure of C(lambda). The corresponding transfer matrices are
Zariski dense in the full image of the map from parameters to transfer
matrices. This establishes failure of algebraically generic
uniqueness even with properness, stability, and the spectral-radius bound.
Genericity in a prescribed positive-order coefficient parameterization
and the strictly proper case require separate analysis.

The manuscript uses "transfer matrix" for C(I-G)^(-1)B, F, and H,
following Legat and Hendrickx, and identifies the corresponding
"closed-loop transfer function" in Hansson, Galvao da Mata, and Andersen.

The abstract contains 192 whitespace-separated words. The compiled paper
has 6 US Letter pages, including the diagram and references, with embedded
Type 1 fonts. All main proofs are included in the manuscript.

The final _algebraic LaTeX log has no undefined references, overfull or underfull
boxes, or package warnings. MiKTeX prints an environment maintenance
notice that an administrator has not checked for updates; this does not
prevent compilation.
