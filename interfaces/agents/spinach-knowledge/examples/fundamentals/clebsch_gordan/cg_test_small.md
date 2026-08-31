# examples/fundamentals/clebsch_gordan/cg_test_small.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/clebsch_gordan/cg_test_small.m`
- Signature: `cg_test_small()`
- Total lines: 35

## Purpose

Compares the output of Spinach Clebsch-Gordan function with the arbitrary precision results returned by Mathematica.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 8-9: Load the test file; implemented by `load('cg_test_table_small.mat','cg_test_table')`.
- Lines 11-12: Loop over the lines of the table; implemented by `for n=1:size(cg_test_table,1)`.
- Lines 14-15: Compute the CG coefficient; implemented by `tic`.
- Lines 21-22: Compare with the Mathematica result; implemented by `difference=abs(cg_test_table(n,7)-cg)`.

### Control flow inferred from the code

- Line 12: `for` loop over `n=1:size(cg_test_table,1)`.
- Line 23: conditional branch on `abs(difference)<2*eps`.

### Key state/data transformations

- Lines 16-18: computes `cg` using `cg=clebsch_gordan(cg_test_table(n,1),cg_test_table(n,2), cg_test_table(n,3),cg_test_table(n,4), cg_test_table(n,5),cg_test_table(n,6))`.
- Lines 22: computes `difference` using `difference=abs(cg_test_table(n,7)-cg)`.

## Implementation structure

- Compares the output of Spinach Clebsch-Gordan function with the
- arbitrary precision results returned by Mathematica.
- Load the test file
- Loop over the lines of the table
- Compute the CG coefficient
- Compare with the Mathematica result

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `clebsch_gordan()`, `cg_test_table()`, `num2str()`.
