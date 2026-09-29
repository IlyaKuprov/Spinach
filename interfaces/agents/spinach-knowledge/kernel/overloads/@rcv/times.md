# kernel/overloads/@rcv/times.m

- Signature: `C=times(A,B)`

## Purpose

Scale one RCV sparse matrix by a numeric scalar, with the scalar in either operand position.

## RCV representation

Here `rcv` means row-column-value: a sparse matrix stored as parallel `row`, `col`, and `val` vectors, with explicit `numRows` and `numCols` dimensions. The class declares the coordinate and dimension vectors as `int64` and values as `double`. This is coordinate-list matrix storage.

Coordinates identify MATLAB matrix row and column positions. The `rcv` class folder defines no custom `subsref` overload; for ordinary element indexing, first convert with `sparse(A)`, then index the MATLAB sparse matrix. That conversion calls `sparse(A.row,A.col,A.val,A.numRows,A.numCols)`. Repeated row-column coordinates can remain as separate stored triplets; MATLAB's sparse constructor combines repeated coordinates by adding their values.

## Inputs

- Exactly one of `A` and `B` is an `rcv` matrix; the other must be a numeric scalar.

## Output

- `C` - the scaled `rcv` matrix.

## Implementation

After validating the operand types and scalar size, the overload multiplies each stored `val` by the scalar and returns the RCV operand. The coordinate vectors `row` and `col` and the explicit dimensions `numRows` and `numCols` are unchanged. Repeated coordinates are scaled independently in storage; later conversion to MATLAB sparse form sums duplicate coordinates. This is scalar scaling, not elementwise multiplication of two matrices or matrix multiplication. For the distinct `mtimes` operation with two RCV operands, the implementation checks `A.numCols == B.numRows`, converts both operands to MATLAB sparse matrices, and returns their MATLAB sparse product, with shape `A.numRows`-by-`B.numCols`.

## Sources

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@rcv/times.m)
- [RCV `mtimes` source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@rcv/mtimes.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=rcv/times.m)
