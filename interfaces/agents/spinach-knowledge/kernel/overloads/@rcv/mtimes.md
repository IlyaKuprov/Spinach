# kernel/overloads/@rcv/mtimes.m

- Signature: `C=mtimes(A,B)`

## Operation and result form

An RCV object stores column-vector coordinates `row` and `col` as `int64`, values in `val` as `double`, dimensions `numRows` and `numCols` as `int64`, and an `isGPU` flag. With an RCV operand and a numeric scalar, the overload scales the RCV `val` vector by that scalar and returns the scaled operand in RCV form. The scalar may be on either side; this is scalar scaling, not array broadcasting. The scalar path does not conjugate values.

For matrix multiplication, the accepted combinations are RCV-by-RCV, RCV-by-MATLAB-sparse, and MATLAB-sparse-by-RCV. Each branch checks the inner dimensions, then converts RCV operand(s) with `sparse` and performs ordinary MATLAB matrix multiplication. The result is a MATLAB sparse matrix with size `size(A,1)-by-size(B,2)`; this path eagerly materialises that sparse product rather than returning a lazy operator. Multiplication uses `*`, not a conjugate-transpose operation.

## Input checks

At least one operand must be RCV. When paired with RCV, the other operand must be RCV, MATLAB sparse, or a numeric scalar; dense nonscalar numeric matrices are not accepted. Scalars are recognised with `isnumeric` and `isscalar`. Matrix cases error when the inner dimensions do not agree. No general implicit expansion or broadcasting is implemented.

## Inputs and outputs

- `A`, `B`: operands in the combinations described above.
- `C`: RCV for scalar scaling; MATLAB sparse for matrix multiplication.

## Source and Wiki

- [Source: kernel/overloads/@rcv/mtimes.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@rcv/mtimes.m)
- [Spinach Wiki: rcv/mtimes.m](https://spindynamics.org/wiki/index.php?title=rcv/mtimes.m)
