# kernel/overloads/@rcv/spy.m

- Signature: `spy(A)`

## Purpose

Plot the nonzero pattern of an RCV sparse matrix.

## RCV representation

Here `rcv` means row-column-value: a sparse matrix stored as parallel `row`, `col`, and `val` vectors, with explicit `numRows` and `numCols` dimensions. The class declares the coordinate and dimension vectors as `int64` and values as `double`. This is coordinate-list matrix storage.

Coordinates identify MATLAB matrix row and column positions. The `rcv` class folder defines no custom `subsref` overload; for ordinary element indexing, first convert with `sparse(A)`, then index the MATLAB sparse matrix. That conversion calls `sparse(A.row,A.col,A.val,A.numRows,A.numCols)`. Repeated row-column coordinates can remain as separate stored triplets; MATLAB's sparse constructor combines repeated coordinates by adding their values.

## Inputs

- `A` - an `rcv` sparse matrix.

## Output

A MATLAB sparsity plot; no matrix is returned.

## Implementation

The overload validates that `A` is an `rcv` object, converts it to MATLAB sparse form with `sparse(A)`, and delegates plotting to MATLAB's `spy`. Thus the plot reflects the matrix after MATLAB has combined any duplicate coordinates, rather than displaying a marker for each raw triplet. Matrix dimensions are retained in the sparse conversion.

## Sources

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@rcv/spy.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=rcv/spy.m)
