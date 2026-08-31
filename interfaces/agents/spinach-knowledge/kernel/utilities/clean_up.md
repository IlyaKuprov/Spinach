# kernel/utilities/clean_up.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/clean_up.m`
- Signature: `A=clean_up(spin_system,A,nonzero_tol)`
- Total lines: 98

## Purpose

Array clean-up utility. Drops non-zero elements with magnitude below the user-specified tolerance and converts between sparse and full storage de- pending on the density of non-zeroes in the array. Syntax: A=clean_up(spin_system,A,nonzero_tol)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 24-25: Skip opium objects; implemented by `if isa(A,'opium'), return; end`.
- Lines 27-28: Skip if disabled; implemented by `if (nonzero_tol==0)||isnan(nonzero_tol), return; end`.
- Lines 30-31: Process cells recursively; implemented by `if iscell(A)`.
- Lines 38-39: Process polyadics recursively; implemented by `if isa(A,'polyadic')`.
- Lines 54-55: Check consistency; implemented by `grumble(A,nonzero_tol)`.
- Lines 57-58: Check if clean-up is allowed; implemented by `if ~ismember('clean-up',spin_system.sys.disable)`.
- Lines 60-61: Use memory-hungry generic method; implemented by `A=nonzero_tol*round((1/nonzero_tol)*A)`.
- Lines 63-65: A small non-zero matrix should always be full; implemented by `if issparse(A)&&(nnz(A)>0)&& any(size(A)<spin_system.tols.small_matrix), A=full(A); end`.
- Lines 67-68: A big sparse matrix with too many non-zeros should be full; implemented by `if issparse(A)&&(nnz(A)/numel(A)>spin_system.tols.dense_matrix), A=full(A); end`.
- Lines 70-72: A big full matrix with too few non-zeros should be sparse; implemented by `if (~issparse(A))&&(nnz(A)/numel(A)<spin_system.tols.dense_matrix)&& (all(size(A)>spin_system.tols.small_matrix))`.

### Control flow inferred from the code

- Line 25: conditional branch on `isa(A,'opium'), return; end`.
- Line 28: conditional branch on `(nonzero_tol==0)||isnan(nonzero_tol), return; end`.
- Line 31: conditional branch on `iscell(A)`.
- Line 32: `for` loop over `n=1:numel(A)`.
- Line 39: conditional branch on `isa(A,'polyadic')`.
- Line 40: `for` loop over `n=1:numel(A.prefix)`.
- Line 43: `for` loop over `n=1:numel(A.suffix)`.
- Line 46: `for` loop over `n=1:numel(A.cores)`.
- Line 47: `for` loop over `k=1:numel(A.cores{n})`.
- Line 58: conditional branch on `~ismember('clean-up',spin_system.sys.disable)`.
- Line 64: conditional branch on `issparse(A)&&(nnz(A)>0)&&`.
- Line 68: conditional branch on `issparse(A)&&(nnz(A)/numel(A)>spin_system.tols.dense_matrix), A=full(A); end`.
- Line 71: conditional branch on `(~issparse(A))&&(nnz(A)/numel(A)<spin_system.tols.dense_matrix)&&`.

### Key state/data transformations

- Lines 33: computes `A{n}` using `A{n}=clean_up(spin_system,A{n},nonzero_tol)`.
- Lines 41: computes `A.prefix{n}` using `A.prefix{n}=clean_up(spin_system,A.prefix{n},nonzero_tol)`.
- Lines 44: computes `A.suffix{n}` using `A.suffix{n}=clean_up(spin_system,A.suffix{n},nonzero_tol)`.
- Lines 48: computes `A.cores{n}{k}` using `A.cores{n}{k}=clean_up(spin_system,A.cores{n}{k},nonzero_tol)`.
- Lines 61: computes `A` using `A=nonzero_tol*round((1/nonzero_tol)*A)`.
- Lines 65: computes `any(size(A)<spin_system.tols.small_matrix), A` using `any(size(A)<spin_system.tols.small_matrix), A=full(A); end`.

### Local helper functions

- Line 81: `grumble()` — `function grumble(A,nonzero_tol)`.
  - Representative operation: `if (~isnumeric(A))`.
  - Representative operation: `error('A must be numeric.')`.

## Parameters / inputs

- A -a numerical array or a cell array thereof
- nonzero_tol -nonzero tolerance

## Outputs

- A -cleaned-up array

## Implementation structure

- Array clean-up utility. Drops non-zero elements with magnitude below the
- user-specified tolerance and converts between sparse and full storage de-
- pending on the density of non-zeroes in the array. Syntax:
- A=clean_up(spin_system,A,nonzero_tol)
- A -a numerical array or a cell array thereof
- nonzero_tol -nonzero tolerance
- A -cleaned-up array
- Skip opium objects
- Skip if disabled
- Process cells recursively
- Process polyadics recursively
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `isnan()`, `iscell()`, `grumble()`, `ismember()`, `issparse()`, `nnz()`, `any()`, `all()`.
