# kernel/utilities/sp_block_diag.m

- Signature: `S=sp_block_diag(varargin)`

## Purpose

Sparse block diagonal matrix from a stack of matrix blocks.

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

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
