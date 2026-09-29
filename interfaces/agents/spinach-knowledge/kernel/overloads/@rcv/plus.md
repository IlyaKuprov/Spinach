# kernel/overloads/@rcv/plus.m

- Signature: `C=plus(A,B)`

## Operation and storage

An RCV object stores column-vector coordinates `row` and `col` as `int64`, values in `val` as `double`, dimensions `numRows` and `numCols` as `int64`, and an `isGPU` flag. For two RCV matrices of equal size, the method implements `A+B` by concatenating their `row`, `col`, and `val` vectors; it does not merge equal coordinates in this method. If either operand is marked GPU-resident, both operands are converted to GPU arrays before concatenation. The output remains RCV and has the common input dimensions. This is an eager construction of the coordinate/value arrays, not a lazy sum operator.

An RCV matrix may also be added to a MATLAB sparse matrix of the same dimensions. The overload checks the sizes, converts the sparse operand to RCV, and recursively uses the RCV-plus-RCV branch. Numeric scalar addition is explicitly rejected: a scalar would make the represented matrix non-sparse. No scalar broadcasting is provided. The implementation adds values directly and does not conjugate them.

## Input checks

At least one operand must be RCV. The other operand must be RCV, MATLAB sparse, or a numeric scalar; the scalar cases then raise the explicit scalar-addition error. Two RCV operands and mixed RCV/sparse operands must have matching row and column dimensions. Other combinations fail the consistency check.

## Inputs and output

- `A`, `B`: equal-sized RCV matrices, or an RCV and a same-sized MATLAB sparse matrix.
- `C`: RCV sum with the common input dimensions.

## Source and Wiki

- [Source: kernel/overloads/@rcv/plus.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@rcv/plus.m)
- [Spinach Wiki: rcv/plus.m](https://spindynamics.org/wiki/index.php?title=rcv/plus.m)
