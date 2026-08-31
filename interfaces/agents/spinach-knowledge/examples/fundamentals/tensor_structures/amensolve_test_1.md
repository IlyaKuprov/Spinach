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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Initialise the random number generator; implemented by `rng(1)`.
- Lines 17-27: Build the main accuracy test cases; implemented by `cases={ struct('name','small_exact','dims',[5 4 3 2],'nterms_A',3,'nterms_y',2,'diag_shift',5e-2,'tol',1e-10, 'opts',struct('nswp',80,'init_guess_rank',2,'enrichment_ran…`.
- Lines 29-30: Run the main dense-reference tests; implemented by `for n=1:numel(cases)`.
- Lines 32-33: Pull out the current case; implemented by `case_data=cases{n}`.
- Lines 35-36: Build the operator, right-hand side, and dense references; implemented by `[A_tt,y_tt,A_dense,y_dense,x_dense,norm_x]=build_spd_case(case_data.dims,case_data.nterms_A,case_data.nterms_y,case_data.diag_shift,n)`.
- Lines 38-39: Run the AMEn solve; implemented by `x_tt=amensolve(A_tt,y_tt,case_data.tol,case_data.opts)`.
- Lines 42-43: Compute the achieved relative solution error; implemented by `rel_err=norm(x_amen-x_dense,2)/max(norm_x,eps)`.
- Lines 45-46: Check the residual against the dense system; implemented by `rel_res=norm(A_dense*x_amen-y_dense,2)/max(norm(y_dense,2),eps)`.
- Lines 48-49: Practical acceptance thresholds for an iterative approximate solver; implemented by `err_lim=max(case_data.err_factor*case_data.tol,1e-10)`.
- Lines 52-53: Check the dense-reference agreement; implemented by `if rel_err>err_lim`.
- Lines 58-59: Check the residual contract; implemented by `if rel_res>res_lim`.
- Lines 64-65: Check that the result is a single tensor train; implemented by `if x_tt.ntrains~=1`.
- Lines 70-71: Check that the physical dimensions are preserved; implemented by `if any(x_tt.sizes(:,2)~=1) || any(x_tt.sizes(:,1)~=A_tt.sizes(:,2))`.
- Lines 76-77: Check that no non-finite numbers appeared; implemented by `if ~all(isfinite(x_amen(:)))`.
- Lines 82-84: Report progress; implemented by `fprintf('amensolve_test_1: %s passed, dim=%d, relerr=%3.3e, relres=%3.3e, maxrank=%d\n', case_data.name,numel(x_dense),rel_err,rel_res,max(x_tt.ranks))`.
- Lines 88-89: Reproducibility check on a medium problem; implemented by `[A_rep,y_rep,~,~,x_rep,norm_rep]=build_spd_case([10 10 10],3,2,8e-2,99)`.
- Lines 102-103: Zero-enrichment regression on a dense-reference system; implemented by `[A_zero,y_zero,A_zero_dense,y_zero_dense,x_zero_dense,norm_zero]=build_spd_case([8 10 10],3,2,8e-2,77)`.
- Lines 114-116: Nonsymmetric smoke test: this path is not the main contract of amensolve, but it should still return a finite answer on a modest dense-reference case.; implemented by `[A_ns_tt,y_ns_tt,A_ns_dense,y_ns_dense,~,~]=build_nonsym_case([6 6 5],4,4,123)`.

### Control flow inferred from the code

- Line 30: `for` loop over `n=1:numel(cases)`.
- Line 53: conditional branch on `rel_err>err_lim`.
- Line 59: conditional branch on `rel_res>res_lim`.
- Line 65: conditional branch on `x_tt.ntrains~=1`.
- Line 71: conditional branch on `any(x_tt.sizes(:,2)~=1) || any(x_tt.sizes(:,1)~=A_tt.sizes(:,2))`.
- Line 77: conditional branch on `~all(isfinite(x_amen(:)))`.
- Line 97: conditional branch on `abs(err_one-err_two)>1e-11`.
- Line 108: conditional branch on `(~isfinite(err_zero))||(~isfinite(res_zero))`.
- Line 120: conditional branch on `(~all(isfinite(x_ns_full(:))))||(~isfinite(res_ns))`.

### Key state/data transformations

- Lines 18-27: computes `cases` using `cases={ struct('name','small_exact','dims',[5 4 3 2],'nterms_A',3,'nterms_y',2,'diag_shift',5e-2,'tol',1e-10, 'opts',struct('nswp',80,'init_guess_rank',2,'enrichment_ran…`.
- Lines 33: computes `case_data` using `case_data=cases{n}`.
- Lines 36: computes `[A_tt,y_tt,A_dense,y_dense,x_dense,norm_x]` using `[A_tt,y_tt,A_dense,y_dense,x_dense,norm_x]=build_spd_case(case_data.dims,case_data.nterms_A,case_data.nterms_y,case_data.diag_shift,n)`.
- Lines 39: computes `x_tt` using `x_tt=amensolve(A_tt,y_tt,case_data.tol,case_data.opts)`.
- Lines 40: computes `x_amen` using `x_amen=full(x_tt)`.
- Lines 43: computes `rel_err` using `rel_err=norm(x_amen-x_dense,2)/max(norm_x,eps)`.
- Lines 46: computes `rel_res` using `rel_res=norm(A_dense*x_amen-y_dense,2)/max(norm(y_dense,2),eps)`.
- Lines 49: computes `err_lim` using `err_lim=max(case_data.err_factor*case_data.tol,1e-10)`.
- Lines 50: computes `res_lim` using `res_lim=max(case_data.res_factor*case_data.tol,1e-10)`.
- Lines 83-84: computes `fprintf('amensolve_test_1: %s passed, dim` using `fprintf('amensolve_test_1: %s passed, dim=%d, relerr=%3.3e, relres=%3.3e, maxrank=%d\n', case_data.name,numel(x_dense),rel_err,rel_res,max(x_tt.ranks))`.
- Lines 89: computes `[A_rep,y_rep,~,~,x_rep,norm_rep]` using `[A_rep,y_rep,~,~,x_rep,norm_rep]=build_spd_case([10 10 10],3,2,8e-2,99)`.
- Lines 90: computes `opts_rep` using `opts_rep=struct('nswp',120,'init_guess_rank',2,'enrichment_rank',6,'rmax',32,'max_full_size',500,'local_iters',220,'verb',0)`.
- Lines 92: computes `x_one` using `x_one=amensolve(A_rep,y_rep,1e-8,opts_rep)`.
- Lines 94: computes `x_two` using `x_two=amensolve(A_rep,y_rep,1e-8,opts_rep)`.
- Lines 95: computes `err_one` using `err_one=norm(full(x_one)-x_rep,2)/max(norm_rep,eps)`.
- Lines 96: computes `err_two` using `err_two=norm(full(x_two)-x_rep,2)/max(norm_rep,eps)`.
- Lines 100: computes `fprintf('amensolve_test_1: reproducibility passed, relerr` using `fprintf('amensolve_test_1: reproducibility passed, relerr=%3.3e\n',err_one)`.
- Lines 104: computes `x_zero` using `x_zero=amensolve(A_zero,y_zero,1e-8,struct('nswp',140,'init_guess_rank',2,'enrichment_rank',0,'rmax',32,'max_full_size',400,'local_iters',220,'verb',0))`.

### Local helper functions

- Line 131: `build_spd_case()` — `function [A_tt,y_tt,A_dense,y_dense,x_dense,norm_x]=build_spd_case(dims,nterms_A,nterms_y,diag_shift,seed)`. Reset the random number generator for this case
  - Representative operation: `rng(seed)`.
  - Representative operation: `B_terms=build_rectangular_terms(dims,nterms_A)`.
- Line 156: `build_nonsym_case()` — `function [A_tt,y_tt,A_dense,y_dense,x_dense,norm_x]=build_nonsym_case(dims,nterms_A,nterms_y,seed)`. Reset the random number generator for this case
  - Representative operation: `rng(seed)`.
  - Representative operation: `A_terms=build_rectangular_terms(dims,nterms_A)`.
- Line 181: `build_rectangular_terms()` — `function terms=build_rectangular_terms(dims,nterms)`. Read the number of cores
  - Representative operation: `ncores=numel(dims)`.
  - Representative operation: `terms=cell(ncores,nterms)`.
- Line 206: `build_vector_terms()` — `function terms=build_vector_terms(dims,nterms)`. Read the number of cores
  - Representative operation: `ncores=numel(dims)`.
  - Representative operation: `terms=cell(ncores,nterms)`.

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
