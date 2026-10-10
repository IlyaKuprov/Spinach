# kernel/overloads/@rcv/horzcat.m

[GitHub source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@rcv/horzcat.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=rcv/horzcat.m)

- Signature: `A=horzcat(varargin)` (the arguments are used in their supplied order)

## Purpose

Horizontally concatenates one or more RCV matrices as consecutive column blocks.

## Storage and behaviour

RCV stores row indices, column indices, and values in parallel arrays, with `numRows` and `numCols` recording the represented shape. The overload requires every argument to be an `rcv` object and all row counts to match; the implementation assumes at least one argument. If any operand has `isGPU=true`, it applies `gpuArray` to every operand, leaving already-GPU operands as they are. In the supplied left-to-right order, each operand's row indices and values are retained, while its column indices are increased by the cumulative widths of all preceding operands. The adjusted row, column, and value arrays are then eagerly concatenated, and `numCols` is set to the sum of the operand widths. The output is an RCV matrix with `numRows` rows and that total number of columns; its data remain in coordinate-array form rather than being materialised as a sparse or dense MATLAB matrix.

Values are concatenated unchanged: this overload does not conjugate them or provide scalar expansion/broadcasting. Non-RCV scalar operands fail the object-type check.

## Input

- One or more RCV sparse matrices, in the order to appear from left to right. Every input must have the same row count.

## Output

- `A` - the RCV matrix formed by the consecutive horizontal blocks.
