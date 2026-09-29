# kernel/overloads/@rcv/transpose.m

- Signature: `A=transpose(A)`

## Purpose

Return the non-conjugating transpose of an RCV sparse matrix.

## RCV representation

Here `rcv` means row-column-value: a sparse matrix stored as parallel `row`, `col`, and `val` vectors, with explicit `numRows` and `numCols` dimensions. The class declares the coordinate and dimension vectors as `int64` and values as `double`. This is coordinate-list matrix storage.

Coordinates identify MATLAB matrix row and column positions. The `rcv` class folder defines no custom `subsref` overload; for ordinary element indexing, first convert with `sparse(A)`, then index the MATLAB sparse matrix. That conversion calls `sparse(A.row,A.col,A.val,A.numRows,A.numCols)`. Repeated row-column coordinates can remain as separate stored triplets; MATLAB's sparse constructor combines repeated coordinates by adding their values.

## Input

- `A` - an `rcv` sparse matrix.

## Output

- `A` - the transposed RCV matrix.

## Implementation

After checking the input class, the overload swaps `row` with `col` and `numRows` with `numCols`. It leaves `val` unchanged, so an entry at `(i,j)` becomes an entry at `(j,i)`, and an `m`-by-`n` matrix becomes `n`-by-`m`. This is `transpose`, not conjugate transpose; use the separate `ctranspose` overload when conjugation is intended. Duplicate coordinates remain duplicated, with coordinates swapped in each triplet.

## Sources

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@rcv/transpose.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=rcv/transpose.m)
