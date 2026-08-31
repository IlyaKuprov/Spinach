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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Initialise the random number generator; implemented by `rng(1)`.
- Lines 20-25: Build the test cases; implemented by `cases={ struct('name','small_exact','dims',[5 4 3 2],'nterms',6,'tol',1e-12,'opts',struct('max_swp',80,'init_guess_rank',2,'enrichment_rank',4,'verb',0)), struct('name',…`.
- Lines 27-28: Run the main numerical tests; implemented by `for n=1:numel(cases)`.
- Lines 30-31: Pull out the current case; implemented by `case_data=cases{n}`.
- Lines 33-34: Build a buffered CP-format tensor train sum and a dense reference; implemented by `[x_ref,dense_ref,norm_ref]=build_case(case_data.dims,case_data.nterms,n)`.
- Lines 36-37: Run the AMEn summation; implemented by `y_ref=amensum(x_ref,case_data.tol,case_data.opts)`.
- Lines 40-41: Compute the achieved relative error; implemented by `rel_err=norm(dense_amen-dense_ref,'fro')/max(norm_ref,eps)`.
- Lines 43-44: Set a practical acceptance threshold for an approximate algorithm; implemented by `rel_lim=max(50*case_data.tol,1e-12)`.
- Lines 46-47: Check the main approximation contract; implemented by `if rel_err>rel_lim`.
- Lines 52-53: Check that the output is a single tensor train; implemented by `if y_ref.ntrains~=1`.
- Lines 58-59: Check that the physical dimensions are preserved; implemented by `if any(any(y_ref.sizes~=x_ref.sizes(:,1:2)))`.
- Lines 64-65: Check that no non-finite numbers appeared; implemented by `if ~all(isfinite(dense_amen(:)))`.
- Lines 70-72: Report progress; implemented by `fprintf('amensum_test_1: %s passed, dim=%d, relerr=%3.3e, maxrank=%d\n', case_data.name,size(dense_ref,1),rel_err,max(y_ref.ranks))`.
- Lines 76-77: Build a reproducibility case for rank-sensitive behaviour; implemented by `[x_rep,dense_rep,norm_rep]=build_case([10 10 10],24,99)`.
- Lines 89-90: Check that disabling enrichment no longer crashes and still returns a finite object; implemented by `[x_smoke,dense_smoke,norm_smoke]=build_case([8 10 10],20,77)`.
- Lines 99-100: Final report; implemented by `fprintf('amensum_test_1: all tests passed.\n')`.

### Control flow inferred from the code

- Line 28: `for` loop over `n=1:numel(cases)`.
- Line 47: conditional branch on `rel_err>rel_lim`.
- Line 53: conditional branch on `y_ref.ntrains~=1`.
- Line 59: conditional branch on `any(any(y_ref.sizes~=x_ref.sizes(:,1:2)))`.
- Line 65: conditional branch on `~all(isfinite(dense_amen(:)))`.
- Line 85: conditional branch on `abs(err_one-err_two)>1e-12`.
- Line 93: conditional branch on `~isfinite(err_smoke)`.

### Key state/data transformations

- Lines 21-25: computes `cases` using `cases={ struct('name','small_exact','dims',[5 4 3 2],'nterms',6,'tol',1e-12,'opts',struct('max_swp',80,'init_guess_rank',2,'enrichment_rank',4,'verb',0)), struct('name',…`.
- Lines 31: computes `case_data` using `case_data=cases{n}`.
- Lines 34: computes `[x_ref,dense_ref,norm_ref]` using `[x_ref,dense_ref,norm_ref]=build_case(case_data.dims,case_data.nterms,n)`.
- Lines 37: computes `y_ref` using `y_ref=amensum(x_ref,case_data.tol,case_data.opts)`.
- Lines 38: computes `dense_amen` using `dense_amen=full(y_ref)`.
- Lines 41: computes `rel_err` using `rel_err=norm(dense_amen-dense_ref,'fro')/max(norm_ref,eps)`.
- Lines 44: computes `rel_lim` using `rel_lim=max(50*case_data.tol,1e-12)`.
- Lines 71-72: computes `fprintf('amensum_test_1: %s passed, dim` using `fprintf('amensum_test_1: %s passed, dim=%d, relerr=%3.3e, maxrank=%d\n', case_data.name,size(dense_ref,1),rel_err,max(y_ref.ranks))`.
- Lines 77: computes `[x_rep,dense_rep,norm_rep]` using `[x_rep,dense_rep,norm_rep]=build_case([10 10 10],24,99)`.
- Lines 78: computes `opts_rep` using `opts_rep=struct('max_swp',100,'init_guess_rank',2,'enrichment_rank',4,'verb',0)`.
- Lines 80: computes `y_one` using `y_one=amensum(x_rep,1e-9,opts_rep)`.
- Lines 82: computes `y_two` using `y_two=amensum(x_rep,1e-9,opts_rep)`.
- Lines 83: computes `err_one` using `err_one=norm(full(y_one)-dense_rep,'fro')/max(norm_rep,eps)`.
- Lines 84: computes `err_two` using `err_two=norm(full(y_two)-dense_rep,'fro')/max(norm_rep,eps)`.
- Lines 90: computes `[x_smoke,dense_smoke,norm_smoke]` using `[x_smoke,dense_smoke,norm_smoke]=build_case([8 10 10],20,77)`.
- Lines 91: computes `y_smoke` using `y_smoke=amensum(x_smoke,1e-8,struct('max_swp',160,'init_guess_rank',2,'enrichment_rank',0,'verb',0))`.
- Lines 92: computes `err_smoke` using `err_smoke=norm(full(y_smoke)-dense_smoke,'fro')/max(norm_smoke,eps)`.
- Lines 96-97: computes `fprintf('amensum_test_1: zero_enrichment smoke passed, dim` using `fprintf('amensum_test_1: zero_enrichment smoke passed, dim=%d, relerr=%3.3e, maxrank=%d\n', size(dense_smoke,1),err_smoke,max(y_smoke.ranks))`.

### Local helper functions

- Line 104: `build_case()` — `function [x_obj,dense_sum,norm_sum]=build_case(dims,nterms,seed)`. Reset the random number generator for this case
  - Representative operation: `rng(seed)`.
  - Representative operation: `ncores=numel(dims)`.

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
