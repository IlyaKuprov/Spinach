# kernel/utilities/kronm_new.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/kronm_new.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/kronm_new.m)

## Purpose

Calculates `(Q{1}(x)Q{2}(x)...(x)Q{n})*M` without opening Kronecker products. This function applies the product of multiple Kronecker factors to a vector or matrix `M` by contracting implicit dimensions individually, avoiding explicit construction of the full Kronecker product.

## Behaviour

- Validates inputs via an internal consistency check (`grumble`):
  - `Q` must be a cell array, otherwise the error `Q must be a cell array.` is raised.
  - Every element of `Q` must be a matrix, otherwise the error `Q must be a cell array of matrices.` is raised.
  - `M` must be numeric, otherwise the error `x must be numeric.` is raised.
- Determines the number of Kronecker terms (`numel(Q)`) and the number of columns of `M` (`size(M,2)`).
- Records row and column dimensions of each factor, taking sizes from the cell in reverse order (`Q{n_mats_in_q-n+1}`), so that `row_dims` and `col_dims` are ordered to match the Kronecker product structure.
- Folds `M` into an implicit multi-dimensional array via `reshape(full(M),[col_dims n_cols_in_m])`, where the leading dimensions correspond to the column dimensions of the Kronecker factors.
- For each factor `n` from 1 to the number of terms, contracts the second dimension of `full(Q{n})` with the corresponding implicit dimension of `M` using `tensorprod(full(Q{n}),M,2,n_mats_in_q)`. Because the contraction index is always `n_mats_in_q` (the last implicit dimension before the column dimensions), each contraction consumes one implicit dimension while the remaining implicit dimensions shift accordingly.
- Unfolds the result back to a matrix with `prod(row_dims)` rows and the original number of columns of `M`.

## Inputs and outputs

**Syntax:** `M=kronm_new(Q,M)`

**Inputs:**

- `Q` — cell array of Kronecker terms (matrices).
- `M` — a vector or a matrix of appropriate dimension.

**Outputs:**

- `M` — a vector or a matrix of appropriate dimension, equal to `(Q{1}(x)Q{2}(x)...(x)Q{n})*M`.

## References

- Spinach wiki page for `kronm.m`: [https://spindynamics.org/wiki/index.php?title=kronm.m](https://spindynamics.org/wiki/index.php?title=kronm.m)
