# kernel/overloads/@rcv/minus.m

- Signature: `A=minus(A,B)`

## Operation

The implementation expresses subtraction in left-to-right order as `A=plus(A,(-1)*B)`: it first scales the right operand by negative one, then passes that result and the left operand to `plus`. It does not swap the operands. The negation is ordinary multiplication by `-1`, with no complex conjugation.

For equal-size RCV matrices, this produces an RCV result with the same `numRows` and `numCols`. The called `plus` method also accepts a MATLAB sparse matrix paired with an RCV matrix, in either order; it checks that their dimensions match, converts the sparse operand to RCV, and returns RCV. Scalar subtraction is not implemented by this path: `plus` rejects adding a numeric scalar to RCV.

## Input checks and storage

The local consistency check errors only when neither operand is an RCV object; it does not itself test dimensions or require both operands to be RCV. When `B` is RCV, negation uses the RCV scalar-by-matrix `mtimes` branch; when `B` is MATLAB sparse, MATLAB performs the scalar multiplication. The subsequent `plus` call enforces equal matrix dimensions for supported matrix operands. In RCV storage, `row` and `col` are `int64` coordinate column vectors, `val` is a `double` value vector, and `numRows`/`numCols` are `int64` dimensions. RCV addition concatenates coordinate/value arrays rather than combining matching coordinates into a single stored entry. The result is an actual RCV object, not a lazy subtraction expression.

## Inputs and output

- `A`: left operand; `B`: right operand. Supported matrix cases include RCV with RCV or one RCV and one MATLAB sparse matrix.
- `A`: result of `A-B` in RCV storage for these supported matrix cases; dimensions are those of the matched operands.

## Source and Wiki

- [Source: kernel/overloads/@rcv/minus.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@rcv/minus.m)
- [Spinach Wiki: rcv/minus.m](https://spindynamics.org/wiki/index.php?title=rcv/minus.m)
