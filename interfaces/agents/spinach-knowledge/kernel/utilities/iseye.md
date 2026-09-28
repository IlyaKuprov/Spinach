# kernel/utilities/iseye.m

- Signature: `verdict=iseye(M)`

## Purpose

Performs the function's computationally affordable test for whether a numeric matrix is the identity matrix.

## Parameters / inputs

- `M` - numeric matrix.

## Outputs

- `verdict` - true or false.

## Numerical / algorithmic content

A nonsquare matrix returns false. For a square matrix, the function first returns false if `M` is not diagonal. Otherwise it draws one random column vector `a` of matching length and returns true exactly when `nnz(M*a-a)` is zero; a failed comparison returns false. The source validates that `M` is numeric.

## Source

[Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=iseye.m)
