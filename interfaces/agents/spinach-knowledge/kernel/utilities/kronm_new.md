# kernel/utilities/kronm_new.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/kronm_new.m`
- Signature: `M=kronm_new(Q,M)`
- Total lines: 71

## Purpose

Calculates (Q{1}(x)Q{2}(x)...(x)Q{n})*M without opening Kronecker products. Syntax: M=kronm(Q,M)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 22-23: Check consistency; implemented by `grumble(Q,M)`.
- Lines 25-26: Dimension statistics; implemented by `n_mats_in_q=numel(Q)`.
- Lines 29-30: Row and column counts in Q; implemented by `row_dims=zeros(1,n_mats_in_q)`.
- Lines 36-37: Fold up implicit dimensions of M; implemented by `M=reshape(full(M),[col_dims n_cols_in_m])`.
- Lines 39-40: Run the products; implemented by `for n=1:n_mats_in_q`.
- Lines 42-43: Contract each implicit dimension; implemented by `M=tensorprod(full(Q{n}),M,2,n_mats_in_q)`.
- Lines 47-48: Unfold implicit dimensions of M; implemented by `M=reshape(M,[prod(row_dims) n_cols_in_m])`.

### Control flow inferred from the code

- Line 32: `for` loop over `n=1:n_mats_in_q`.
- Line 40: `for` loop over `n=1:n_mats_in_q`.

### Key state/data transformations

- Lines 26: computes `n_mats_in_q` using `n_mats_in_q=numel(Q)`.
- Lines 27: computes `n_cols_in_m` using `n_cols_in_m=size(M,2)`.
- Lines 30: computes `row_dims` using `row_dims=zeros(1,n_mats_in_q)`.
- Lines 31: computes `col_dims` using `col_dims=zeros(1,n_mats_in_q)`.
- Lines 33: computes `[row_dims(n),col_dims(n)]` using `[row_dims(n),col_dims(n)]=size(Q{n_mats_in_q-n+1})`.
- Lines 37: computes `M` using `M=reshape(full(M),[col_dims n_cols_in_m])`.

### Local helper functions

- Line 53: `grumble()` — `function grumble(Q,x)`.
  - Representative operation: `if (~iscell(Q))`.
  - Representative operation: `error('Q must be a cell array.')`.

## Parameters / inputs

- Q -cell array of Kronecker terms
- M -a vector or a matrix of appropriate dimension
- Output:
- M -a vector or a matrix of appropriate dimension

## Implementation structure

- Calculates (Q{1}(x)Q{2}(x)...(x)Q{n})*M without opening
- Kronecker products. Syntax:
- M=kronm(Q,M)
- Q - cell array of Kronecker terms
- M - a vector or a matrix of appropriate dimension
- Output:
- Check consistency
- Dimension statistics
- Row and column counts in Q
- Fold up implicit dimensions of M
- Run the products
- Contract each implicit dimension

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `row_dims()`, `col_dims()`, `tensorprod()`, `iscell()`, `ismatrix()`.
