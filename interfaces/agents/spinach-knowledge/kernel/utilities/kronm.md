# kernel/utilities/kronm.m

## Purpose

`kronm.m` calculates `(Q{1}(x)Q{2}(x)...(x)Q{n})*x` — the action of a Kronecker product of matrices on a vector or matrix — without explicitly opening (forming) the Kronecker products.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/kronm.m>

## Behaviour

- Syntax: `x=kronm(Q,x)`.
- A consistency check (`grumble`) is run first: `Q` must be a cell array, every element of `Q` must be a matrix, and `x` must be numeric; otherwise errors are thrown (`'Q must be a cell array.'`, `'Q must be a cell array of matrices.'`, `'x must be numeric.'`).
- The number of matrices in `Q` is `nmats=numel(Q)`; the number of columns of `x` is `ncols=size(x,2)`.
- Row and column counts of each factor are collected in reverse order: for `n=1:nmats`, `[row_dims(n),col_dims(n)]=size(Q{nmats-n+1})`.
- A dimension map for `x` is built as `x_dims=[col_dims,ncols]`, and `x` is reshaped (after `full`) into that map.
- The products run over `n=1:nmats`:
  - Shortcut for `opium` objects: if `isa(Q{nmats-n+1},'opium')` and its `coeff` is not equal to 1, the step is `x=Q{nmats-n+1}.coeff*x` and the loop continues.
  - For `n==1` (the leading dimension), no permutation is needed: `x` is reshaped to `[x_dims(1) prod(x_dims)/x_dims(1)]`, multiplied as `x=Q{nmats}*x`, the dimension map is updated with `x_dims(1)=row_dims(1)`, and `x` is reshaped back with `full`.
  - Otherwise, `permute` is used: the `n`-th dimension is brought forward via `dims=[n,setdiff(1:numel(x_dims),n)]`, `x` is reshaped to `[col_dims(n),numel(x)/col_dims(n)]`, multiplied as `x=Q{nmats-n+1}*x`, the dimension map is updated with `x_dims(n)=row_dims(n)`, `x` is reshaped to `[row_dims(n),x_dims(dims(2:end))]`, and `ipermute` restores the dimension order.
- Finally, `x` is reshaped to `[prod(row_dims),ncols]` for output.

## Inputs and outputs

**Inputs**

- `Q` — cell array of Kronecker terms (each element a matrix; `opium` objects are handled via the coefficient shortcut).
- `x` — a vector or a matrix of appropriate dimension; must be numeric.

**Output**

- `x` — a vector or a matrix of appropriate dimension, the result of the Kronecker-product action.

## References

- Spinach Wiki page for `kronm.m`: <https://spindynamics.org/wiki/index.php?title=kronm.m>
- Source file: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/kronm.m>
