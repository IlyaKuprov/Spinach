# examples/fundamentals/tensor_structures/amensolve_test_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/tensor_structures/amensolve_test_1.m`
- Signature: `amensolve_test_1()`
- Total lines: 226

## Purpose

Detailed unit test for ttclass/amensolve against dense references. The test builds structured positive-definite tensor-train linear systems, solves them with AMEn, and compares the result against dense direct solves and dense residuals. The cases include exact small systems, a dense-reference case of dimension 2000, a zero-enrichment regression, and a nonsymmetric finite-output smoke.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- The file also defines local helper function(s): `build_spd_case()`, `build_nonsym_case()`, `build_rectangular_terms()`, `build_vector_terms()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Implementation structure

- Detailed unit test for ttclass/amensolve against dense references.
- The test builds structured positive-definite tensor-train linear systems,
- solves them with AMEn, and compares the result against dense direct solves
- and dense residuals. The cases include exact small systems, a dense-reference
- case of dimension 2000, a zero-enrichment regression, and a nonsymmetric
- finite-output smoke.
- Initialise the random number generator
- Build the main accuracy test cases
- Run the main dense-reference tests
- Pull out the current case
- Build the operator, right-hand side, and dense references
- Run the AMEn solve

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `rng()`, `build_spd_case()`, `amensolve()`, `any()`, `all()`, `x_amen()`, `build_nonsym_case()`, `x_ns_full()`, `build_rectangular_terms()`, `build_vector_terms()`, `ttclass()`, `shrink()`, `unit_like()`, `dims()`.
