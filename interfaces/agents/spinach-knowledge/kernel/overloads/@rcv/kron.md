# kernel/overloads/@rcv/kron.m

- Signature: `C=kron(A,B)`

## Operation

This overload forms the Kronecker product of two RCV matrices. If `A` has dimensions `m-by-n` and `B` has dimensions `p-by-q`, the result has dimensions `(m*p)-by-(n*q)`. Each stored pair `A(i,j)`, `B(k,l)` contributes `A(i,j)*B(k,l)` at row `(i-1)*p+k`, column `(j-1)*q+l`.

## RCV storage and execution

An RCV object stores column-vector coordinates `row` and `col` as `int64`, values in `val` as `double`, dimensions `numRows` and `numCols` as `int64`, and an `isGPU` flag; when GPU-resident, these arrays are GPU arrays. This method requires both operands to be RCV objects. If either is marked GPU-resident, it converts both to GPU arrays before forming `ndgrid(1:length(A.val),1:length(B.val))`, mapping every stored-entry pair, multiplying the paired values, and constructing the output RCV object. The coordinate/value arrays for all pairs are built eagerly; the result is not a lazy Kronecker operator. The output GPU flag is the logical OR of the input flags.

The value operation is ordinary elementwise multiplication; it does not conjugate either operand. There is no scalar or broadcast path in this overload.

## Inputs and output

- `A` and `B`: RCV matrices; the consistency check rejects a non-RCV operand.
- `C`: RCV matrix of size `(A.numRows*B.numRows)-by-(A.numCols*B.numCols)`.

## Source and Wiki

- [Source: kernel/overloads/@rcv/kron.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@rcv/kron.m)
- [Spinach Wiki: rcv/kron.m](https://spindynamics.org/wiki/index.php?title=rcv/kron.m)
