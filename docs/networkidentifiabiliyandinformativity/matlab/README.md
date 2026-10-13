# Reproducing the computational examples

The five MATLAB files in this directory reproduce the two examples in the manuscript's Computational Examples section. They require MATLAB R2025a (25.1) and Symbolic Math Toolbox R2025a (25.1). The scripts use symbolic parameters for the unknown transfer functions; no numerical data or fitted model is involved.

From this directory, run:

```matlab
test_rational_dimension_cycle_5_c
test_rational_dimension_cycle_5_c_subnetwork
```

The first script should report image dimension 11 and generic full row rank (`case 3`, with exceptional rank loss at `c1 = 0`, corresponding to `R^1 = 0`). The second should report image dimension 4 with either three or two active measured inputs, and generic full row rank with five, three, or two active measured inputs (`case 3` in each case). In these rank tests, `case 3` means full row rank holds generically but fails at exceptional parameter values.

The first image-dimension calculation took 630.54 seconds in one run on 7 October 2026 with an AMD Ryzen 7 250 processor (8 cores, 16 logical processors), 31.3 GiB of physical memory, and the MATLAB and toolbox versions above. This timing covers only `dimMatrixRationalImage`, as measured by the script's `tic`/`toc`; MATLAB startup and the rank test are separate. Runtime depends on the machine and MATLAB release.

The two `test_rational_dimension_*.m` files are entry points. `dimMatrixRationalImage.m` computes image dimension by elimination; `findMaxIndependentSet.m` computes dimension from the leading monomials of the elimination basis; `groebner_rank_test_M4.m` tests generic full row rank on the admissible domain.
