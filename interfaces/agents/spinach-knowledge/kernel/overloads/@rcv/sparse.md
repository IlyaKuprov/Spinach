# kernel/overloads/@rcv/sparse.m

- Signature: `A=sparse(A)`
- Source: [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@rcv/sparse.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=rcv/sparse.m)

## Purpose

Converts an RCV coordinate-storage object to an ordinary MATLAB sparse matrix with the same dimensions.

## Behaviour

The function requires an RCV input. If there are no stored entries, it creates a sparse matrix with `numRows` rows and `numCols` columns. Otherwise it calls MATLAB's `sparse(row,col,val,numRows,numCols)` constructor. Consequently, repeated RCV coordinates contribute their values additively at that matrix position in the converted matrix.

The output is a MATLAB sparse matrix, not an RCV object. Its shape comes from the stored dimensions, including when the coordinate arrays are empty.
