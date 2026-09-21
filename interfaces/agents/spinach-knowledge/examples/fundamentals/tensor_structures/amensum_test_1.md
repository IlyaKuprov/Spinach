# examples/fundamentals/tensor_structures/amensum_test_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/tensor_structures/amensum_test_1.m`
- Signature: `amensum_test_1()`
- Total lines: 150

## Purpose

Detailed unit test for ttclass/amensum against dense references. The test uses buffered rank-one tensor trains, compares the AMEn summation result against the exact dense sum, and checks both relative Frobenius error and internal consistency properties. Note: the underlying AMEn summation is approximate, and the paper motivating the method focuses on enrichment-assisted updates. The strict accuracy checks below there

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- The file also defines local helper function(s): `build_case()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Implementation structure

- Detailed unit test for ttclass/amensum against dense references.
- The test uses buffered rank-one tensor trains, compares the AMEn
- summation result against the exact dense sum, and checks both
- relative Frobenius error and internal consistency properties.
- Note: the underlying AMEn summation is approximate, and the paper
- motivating the method focuses on enrichment-assisted updates. The
- strict accuracy checks below therefore target enriched runs, while
- zero-enrichment is kept as a finite-output regression smoke test.
- Initialise the random number generator
- Build the test cases
- Run the main numerical tests
- Pull out the current case

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `rng()`, `build_case()`, `amensum()`, `any()`, `all()`, `dense_amen()`, `coeff()`, `dims()`, `ttclass()`.
