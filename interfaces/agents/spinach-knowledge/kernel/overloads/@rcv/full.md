# kernel/overloads/@rcv/full.m

[GitHub source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@rcv/full.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=rcv/full.m)

- Signature: `A=full(A)`

## Purpose

Materialises the matrix encoded by an RCV sparse matrix as a dense MATLAB matrix at its recorded dimensions.

## Storage and behaviour

RCV stores row indices, column indices, and corresponding values in parallel arrays, with `numRows` and `numCols` retaining the matrix shape. This overload checks that the input is an `rcv` object, then evaluates `full(sparse(A))`: the RCV sparse conversion first constructs a MATLAB sparse matrix at the recorded dimensions, and MATLAB's `full` then eagerly creates the dense result. The output is a full MATLAB matrix of size `numRows`-by-`numCols`; this is materialisation, not a lazy RCV result.

The overload does not conjugate values or implement scalar expansion or broadcasting.

## Input

- `A` - an RCV sparse matrix. The explicit check is object type only.

## Output

- `A` - the corresponding full MATLAB matrix.
