# kernel/overloads/@rcv/ctranspose.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@rcv/ctranspose.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=rcv/ctranspose.m)

- Signature: `A=ctranspose(A)`

## Purpose and representation

Returns the conjugate transpose of an RCV sparse matrix while keeping its coordinate-style storage: row indices in `A.row`, column indices in `A.col`, values in `A.val`, and shape in `A.numRows` and `A.numCols`.

## Operation and dimensions

The implementation first checks that `A` is an `rcv` object. It then swaps `A.row` with `A.col`, swaps `A.numRows` with `A.numCols`, and replaces `A.val` with `conj(A.val)`. Thus each stored coordinate/value pair `(i,j,v)` becomes `(j,i,conj(v))`, and the output shape is `[numCols,numRows]` from the input shape.

Conjugation is elementwise over the stored value array. The implementation does not materialise a full matrix and does not add explicit checks for coordinate/value-array consistency or scalar/broadcast behaviour.
