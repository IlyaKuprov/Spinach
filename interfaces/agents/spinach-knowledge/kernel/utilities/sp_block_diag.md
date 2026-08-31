# kernel/utilities/sp_block_diag.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/sp_block_diag.m`
- Signature: `S=sp_block_diag(varargin)`
- Total lines: 109

## Purpose

Sparse block diagonal matrix from a stack of matrix blocks.

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 29-30: Check consistency; implemented by `if nargin==0`.
- Lines 35-36: Single input: block stack; implemented by `if nargin==1`.
- Lines 38-39: Get the array; implemented by `A=varargin{1}`.
- Lines 41-42: Get block dimensions; implemented by `nrows=size(A,1); ncols=size(A,2)`.
- Lines 50-51: Reshape higher dimensions into a block index; implemented by `A=reshape(A,[nrows ncols nblocks])`.
- Lines 53-54: Get non-zero elements; implemented by `[lin_idx,~,vals]=find(A(:))`.
- Lines 57-59: Assemble the matrix; implemented by `S=sparse(rows+(blocks-1)*nrows,cols+(blocks-1)*ncols,vals, nrows*nblocks,ncols*nblocks)`.
- Lines 63-64: Count block dimensions; implemented by `nrows=zeros(nargin,1); ncols=zeros(nargin,1)`.
- Lines 70-71: Compute offsets; implemented by `row_offsets=[0; cumsum(nrows(1:(end-1)))]`.
- Lines 74-75: Preallocate index cells; implemented by `rows=cell(nargin,1); cols=cell(nargin,1); vals=cell(nargin,1)`.
- Lines 77-78: Loop over blocks; implemented by `for n=1:nargin`.
- Lines 85-87: Assemble the matrix; implemented by `S=sparse(vertcat(rows{:}),vertcat(cols{:}),vertcat(vals{:}), sum(nrows),sum(ncols))`.

### Control flow inferred from the code

- Line 30: conditional branch on `nargin==0`.
- Line 36: conditional branch on `nargin==1`.
- Line 44: conditional branch on `numel(block_dims)<3`.
- Line 65: `for` loop over `n=1:nargin`.
- Line 78: `for` loop over `n=1:nargin`.

### Key state/data transformations

- Lines 39: computes `A` using `A=varargin{1}`.
- Lines 42: computes `nrows` using `nrows=size(A,1); ncols=size(A,2)`.
- Lines 43: computes `block_dims` using `block_dims=size(A)`.
- Lines 45: computes `nblocks` using `nblocks=1`.
- Lines 54: computes `[lin_idx,~,vals]` using `[lin_idx,~,vals]=find(A(:))`.
- Lines 55: computes `[rows,cols,blocks]` using `[rows,cols,blocks]=ind2sub([nrows ncols nblocks],lin_idx)`.
- Lines 58-59: computes `S` using `S=sparse(rows+(blocks-1)*nrows,cols+(blocks-1)*ncols,vals, nrows*nblocks,ncols*nblocks)`.
- Lines 66: computes `nrows(n)` using `nrows(n)=size(varargin{n},1)`.
- Lines 67: computes `ncols(n)` using `ncols(n)=size(varargin{n},2)`.
- Lines 71: computes `row_offsets` using `row_offsets=[0; cumsum(nrows(1:(end-1)))]`.
- Lines 72: computes `col_offsets` using `col_offsets=[0; cumsum(ncols(1:(end-1)))]`.
- Lines 75: computes `rows` using `rows=cell(nargin,1); cols=cell(nargin,1); vals=cell(nargin,1)`.
- Lines 79: computes `[rows{n},cols{n},vals{n}]` using `[rows{n},cols{n},vals{n}]=find(varargin{n})`.
- Lines 80: computes `rows{n}` using `rows{n}=rows{n}(:)+row_offsets(n)`.
- Lines 81: computes `cols{n}` using `cols{n}=cols{n}(:)+col_offsets(n)`.
- Lines 82: computes `vals{n}` using `vals{n}=vals{n}(:)`.

### Local helper functions

- Line 94: `grumble()` — `function grumble(varargin)`. Hopeless? It's not hopeless.
  - Representative operation: `for n=1:nargin`.
  - Representative operation: `if (~isfloat(varargin{n}))||(~isnumeric(varargin{n}))`.

## Syntax

```matlab
S=sp_block_diag(A)
S=sp_block_diag(A,B,C,...)
```

## Parameters / inputs

- A -a floating point array; in the first syntax, the
- first two dimensions are matrix dimensions and the
- remaining dimensions enumerate the blocks
- A,B,C -floating point matrices to be placed on the block
- diagonal in the second syntax

## Outputs

- S -sparse block diagonal matrix
- Notes: this function is a Spinach-local replacement for Matlab's
- spblkdiag function from the Model-Based Calibration toolbox.

## Implementation structure

- Sparse block diagonal matrix from a stack of matrix blocks.
- S=sp_block_diag(A)
- S=sp_block_diag(A,B,C,...)
- A -a floating point array; in the first syntax, the
- first two dimensions are matrix dimensions and the
- remaining dimensions enumerate the blocks
- A,B,C -floating point matrices to be placed on the block
- diagonal in the second syntax
- S -sparse block diagonal matrix
- spblkdiag function from the Model-Based Calibration toolbox.
- Check consistency
- Single input: block stack

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `block_dims()`, `ind2sub()`, `nrows()`, `ncols()`, `cumsum()`, `row_offsets()`, `col_offsets()`, `vertcat()`, `isfloat()`, `ismatrix()`.
