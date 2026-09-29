# kernel/overloads/@ttclass/vec.m

Source: [kernel/overloads/@ttclass/vec.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/vec.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=ttclass/vec.m)

- Signature: `A=vec(A)`

## Purpose

Reshape an ordinary array to a column vector, or collapse the two physical dimensions of every tensor-train core without assembling the full tensor train.

## Behaviour

For non-`ttclass` input, the implementation returns `reshape(A,[numel(A) 1])`. MATLAB's column-major linear ordering is retained, and the output shape is `numel(A)`-by-1.

For a `ttclass`, the method reads the core count and train count from `size(A.cores)`, and obtains the core ranks and two physical dimensions from `ranks(A)` and `sizes(A)`. It visits each train column `n` and each core row `k`, replacing `A.cores{k,n}` with a reshape of dimensions

`[ttm_ranks(k,n), ttm_sizes(k,1)*ttm_sizes(k,2), 1, ttm_ranks(k+1,n)]`.

The core-cell array's core/train ordering is retained. Within each core, MATLAB reshape preserves linear order: physical indices `(i,j)` are packed as `i + (j-1)*ttm_sizes(k,1)`. The left and right bond ranks are unchanged; the two physical dimensions become one dimension and a singleton third dimension. The returned `A` is still a tensor-train object, with the core cells replaced by their reshaped values. This is not the same as column-wise reshaping the fully assembled tensor: the element order can differ by a permutation, as the source warning notes, while matching the tensor-train Kronecker-product ordering.

No rank truncation, rounding, or tolerance test is performed: this operation only reshapes core arrays, leaving the ranks unchanged.

## Input and output

- Input `A`: an ordinary array or a `ttclass` array.
- Output `A`: a column vector for an ordinary array; a `ttclass` result whose individual cores have the shapes described above.
