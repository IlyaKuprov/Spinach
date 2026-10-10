# kernel/utilities/sp_block_diag.m

## Purpose

Builds a sparse block diagonal matrix, either from a multidimensional stack of matrix blocks or from a list of separate matrices. The function is a Spinach-local replacement for MATLAB's `spblkdiag` function from the Model-Based Calibration toolbox.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/sp_block_diag.m>

## Behaviour

Two syntaxes are supported:

- `S=sp_block_diag(A)` — the first two dimensions of `A` are treated as matrix dimensions and all remaining dimensions enumerate the blocks. If `A` has fewer than three dimensions, a single block is used. Higher dimensions are reshaped into a single block index, and the non-zero elements of each block are placed at positions offset by `(block-1)*nrows` in rows and `(block-1)*ncols` in columns, producing an `nrows*nblocks` by `ncols*nblocks` sparse matrix.
- `S=sp_block_diag(A,B,C,...)` — each input matrix is placed on the block diagonal. Row and column offsets are computed as cumulative sums of the preceding blocks' dimensions, and the non-zero elements of all blocks are assembled into a sparse matrix of size `sum(nrows)` by `sum(ncols)`.

Input validation (via the local `grumble` subfunction):

- At least one input array must be supplied; otherwise an error is thrown.
- All inputs must be floating-point numeric arrays (`double` or `single`); otherwise an error is thrown.
- With multiple inputs, every input must be a matrix (2-D); otherwise an error is thrown.

## Inputs and outputs

**Inputs**

- `A` — floating-point array. In the single-input syntax, the first two dimensions are matrix dimensions and the remaining dimensions enumerate the blocks.
- `A,B,C,...` — floating-point matrices to be placed on the block diagonal in the multiple-input syntax.

**Outputs**

- `S` — sparse block diagonal matrix.

## References

- Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=sp_block_diag.m>
- Source file: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/sp_block_diag.m>
