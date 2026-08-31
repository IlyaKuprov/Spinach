# kernel/utilities/kronm.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/kronm.m`
- Signature: `x=kronm(Q,x)`
- Total lines: 111

## Purpose

Calculates (Q{1}(x)Q{2}(x)...(x)Q{n})*x without opening Kronecker products. Syntax: x=kronm(Q,x)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 22-23: Check consistency; implemented by `grumble(Q,x)`.
- Lines 25-26: Number of matrices in Q; implemented by `nmats=numel(Q)`.
- Lines 28-29: Number of columns in x; implemented by `ncols=size(x,2)`.
- Lines 31-32: Row and column counts in Q; implemented by `row_dims=zeros(1,nmats)`.
- Lines 38-39: Dimension map for x; implemented by `x_dims=[col_dims,ncols]`.
- Lines 41-42: Reshape into the map; implemented by `x=reshape(full(x),x_dims)`.
- Lines 44-45: Run the products; implemented by `for n=1:nmats`.
- Lines 47-48: Shortcut for opium objects; implemented by `if isa(Q{nmats-n+1},'opium')&&(Q{nmats-n+1}.coeff~=1)`.
- Lines 52-53: If our dimension is first, do not permute; implemented by `if n==1`.
- Lines 55-56: Unroll other dimensions; implemented by `x=reshape(x,[x_dims(1) prod(x_dims)/x_dims(1)])`.
- Lines 58-59: Run multiplication and update dimension map; implemented by `x=Q{nmats}*x; x_dims(1)=row_dims(1)`.
- Lines 61-62: Roll other dimensions back up; implemented by `x=reshape(full(x),x_dims)`.
- Lines 64-65: Otherwise, permute() is unavoidable; implemented by `else`.
- Lines 67-68: Bring forward n-th dimension; implemented by `dims=1:numel(x_dims); dims(n)=[]`.
- Lines 71-72: Unroll other dimensions; implemented by `x=reshape(x,[col_dims(n),numel(x)/col_dims(n)])`.
- Lines 74-75: Run multiplication and update dimension map; implemented by `x=Q{nmats-n+1}*x; x_dims(n)=row_dims(n)`.
- Lines 77-78: Roll other dimensions back up; implemented by `x=reshape(full(x),[row_dims(n),x_dims(dims(2:end))])`.
- Lines 80-81: Put the current dimension back; implemented by `x=ipermute(x,dims)`.

### Control flow inferred from the code

- Line 34: `for` loop over `n=1:nmats`.
- Line 45: `for` loop over `n=1:nmats`.
- Line 48: conditional branch on `isa(Q{nmats-n+1},'opium')&&(Q{nmats-n+1}.coeff~=1)`.
- Line 53: conditional branch on `n==1`.

### Key state/data transformations

- Lines 26: computes `nmats` using `nmats=numel(Q)`.
- Lines 29: computes `ncols` using `ncols=size(x,2)`.
- Lines 32: computes `row_dims` using `row_dims=zeros(1,nmats)`.
- Lines 33: computes `col_dims` using `col_dims=zeros(1,nmats)`.
- Lines 35: computes `[row_dims(n),col_dims(n)]` using `[row_dims(n),col_dims(n)]=size(Q{nmats-n+1})`.
- Lines 39: computes `x_dims` using `x_dims=[col_dims,ncols]`.
- Lines 42: computes `x` using `x=reshape(full(x),x_dims)`.
- Lines 68: computes `dims` using `dims=1:numel(x_dims); dims(n)=[]`.

### Local helper functions

- Line 93: `grumble()` — `function grumble(Q,x)`.
  - Representative operation: `if (~iscell(Q))`.
  - Representative operation: `error('Q must be a cell array.')`.

## Parameters / inputs

- Q -cell array of Kronecker terms
- x -a vector or a matrix of appropriate dimension
- Output:
- x -a vector or a matrix of appropriate dimension

## Implementation structure

- Calculates (Q{1}(x)Q{2}(x)...(x)Q{n})*x without opening
- Kronecker products. Syntax:
- x=kronm(Q,x)
- Q - cell array of Kronecker terms
- x - a vector or a matrix of appropriate dimension
- Output:
- Check consistency
- Number of matrices in Q
- Number of columns in x
- Row and column counts in Q
- Dimension map for x
- Reshape into the map

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `row_dims()`, `col_dims()`, `x_dims()`, `dims()`, `ipermute()`, `iscell()`, `ismatrix()`.
