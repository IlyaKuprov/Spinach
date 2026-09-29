# kernel/overloads/@rcv/vertcat.m

- Signature: `A=vertcat(varargin)`

## Purpose

Vertically concatenate two or more RCV sparse matrices in operand order.

## RCV representation

Here `rcv` means row-column-value: a sparse matrix stored as parallel `row`, `col`, and `val` vectors, with explicit `numRows` and `numCols` dimensions. The class declares the coordinate and dimension vectors as `int64` and values as `double`. This is coordinate-list matrix storage.

Coordinates identify MATLAB matrix row and column positions. The `rcv` class folder defines no custom `subsref` overload; for ordinary element indexing, first convert with `sparse(A)`, then index the MATLAB sparse matrix. That conversion calls `sparse(A.row,A.col,A.val,A.numRows,A.numCols)`. Repeated row-column coordinates can remain as separate stored triplets; MATLAB's sparse constructor combines repeated coordinates by adding their values.

## Inputs

- One or more `rcv` matrices; all must have the same `numCols`.

## Output

- `A` - the concatenated RCV matrix with row count equal to the sum of the input row counts and the shared column count.

## Implementation

The overload rejects non-RCV inputs and mismatched column counts. If any operand is on the GPU, it converts all operands to GPU arrays before combining them. For each input in order, it adds the cumulative number of preceding rows to that input's `row` coordinates; it copies the `col` and `val` vectors unchanged, then concatenates each vector once. The result's `numRows` is the sum of the inputs' row counts; `numCols` stays the common input column count. Stored duplicates are retained as triplets, including when coordinates from different blocks overlap in column but their row offsets differ.

## Sources

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@rcv/vertcat.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=rcv/vertcat.m)
