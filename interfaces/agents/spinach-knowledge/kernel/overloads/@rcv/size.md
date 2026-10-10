# kernel/overloads/@rcv/size.m

- Signature: `s=size(A,dim)` or `[s,ncols]=size(A)`
- Source: [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@rcv/size.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=rcv/size.m)

## Purpose

Returns the stored matrix dimensions of an RCV sparse-matrix object.

## Behaviour

- With one output and no dimension index, returns `[A.numRows A.numCols]`.
- With one output and `dim=1` or `dim=2`, returns the row or column count, respectively. Any other positive integer dimension returns `1`.
- With two outputs and no dimension index, returns the row count in `s` and the column count in `ncols`.
- A dimension index must be a finite, real, positive integer scalar. Requesting two outputs together with a dimension index raises an error.

Dimensions are carried in the RCV object as `int64` values. The returned dimension answers depend on `numRows` and `numCols`, not the number of stored coordinate entries.
