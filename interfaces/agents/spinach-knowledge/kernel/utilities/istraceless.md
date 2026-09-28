# kernel/utilities/istraceless.m

- Signature: `A=istraceless(M)`

## Purpose

Tests whether a numeric matrix is traceless within a tolerance set by the floating-point precision of its class and the matrix norm.

## Parameters / inputs

- `M` - numeric matrix of any dimension.

## Outputs

- `A` - true when the trace test passes; false otherwise.

## Numerical / algorithmic content

The function sets `precision=eps(class(M))`, computes `norm_m=cheap_norm(M)`, and returns whether `abs(trace(M)) <= precision*norm_m`. It rejects nonnumeric inputs.

## Source

[Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=istraceless.m)
