# kernel/utilities/sp_block_diag.m

- Signature: `S=sp_block_diag(varargin)`

## Purpose

Construct a sparse block diagonal matrix from a stack of matrix blocks or from separate matrices.

## Syntax

```matlab
S=sp_block_diag(A)
S=sp_block_diag(A,B,C,...)
```

## Parameters / inputs

- `A` — a floating-point array. In the single-input syntax, its first two dimensions are matrix dimensions and any remaining dimensions enumerate the blocks.
- `A,B,C,...` — floating-point matrices to place on the block diagonal in the multiple-input syntax.
- At least one input is required. Inputs must be single- or double-precision arrays; with multiple inputs, each must be a matrix.

## Outputs

- `S` — sparse block diagonal matrix.

## Numerical / algorithmic content

- With one input, higher dimensions of `A` are reshaped into a block index. Nonzero entries are placed using row and column offsets based on the block dimensions.
- With multiple inputs, the row and column offsets are cumulative sums of the preceding blocks’ dimensions. Nonzero entries from each matrix are placed at those offsets.

## Notes

This function is a Spinach-local replacement for Matlab's `spblkdiag` function from the Model-Based Calibration toolbox.

<https://spindynamics.org/wiki/index.php?title=sp_block_diag.m>

ilya.kuprov@weizmann.ac.il