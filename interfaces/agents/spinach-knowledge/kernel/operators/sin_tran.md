# kernel/operators/sin_tran.m

- Source: [kernel/operators/sin_tran.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/operators/sin_tran.m)
- Wiki: [sin_tran.m on the Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=sin_tran.m)
- Signature: `A=sin_tran(dim)`

## Purpose

Returns the single-transition matrix units spanning the full space of `dim`-by-`dim` matrices. This is an operator basis, not a generator or a propagator, and the implementation does not assign it a spin or boson model.

## Ordering and construction

The output is a `dim^2`-by-1 cell array. Each entry is a `dim`-by-`dim` sparse complex matrix with exactly one nonzero value, equal to 1, at its assigned row and column. Cell indices enumerate matrix entries along successive anti-diagonals: start at the first column and move upward along each anti-diagonal, then continue to the next one. For dimension 4, the cell index at each matrix location is:

```text
 1   3   6  10
 2   5   9  13
 4   8  12  15
 7  11  14  16
```

Thus each cell contains a matrix unit, and the `dim^2` cells cover all row/column positions exactly once. The implementation computes each position with `lin2kq(dim,n,1)` and constructs the corresponding sparse unit-entry matrix; it declares the result complex and uses `parfor` over the cell indices.

## Input

- `dim`: positive real integer giving the row and column dimension.

## Output

- `A`: a column cell array of `dim^2` sparse complex matrix units in the ordering above.
