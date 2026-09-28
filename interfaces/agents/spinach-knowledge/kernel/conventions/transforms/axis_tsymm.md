# kernel/conventions/transforms/axis_tsymm.m

- Signature: `A=axis_tsymm(T,a)`

## Purpose

Approximately average a real interaction tensor over rotations about a specified axis.

## Physical / mathematical content

For each sampled rotation matrix `R` about axis `a`, the tensor is transformed as `R*T*R'`. The output `A` is the mean of these transformed tensors.

## Numerical / algorithmic content

The function samples 360 rotations, at angles `pi*n/180` for `n=1:360`, and divides their sum by 360.

## Parameters / inputs

- T -3x3 real interaction tensor
- a -3x1 real vector specifying the rotation axis

## Outputs

- A -3x3 real interaction tensor averaged over
- the rotation around the specified axis

## Implementation structure

The function checks that `a` is a real, nonzero 3x1 numeric vector and `T` is a real 3x3 numeric matrix. It initializes a 3x3 zero matrix, accumulates the rotated tensors using `anax2dcm`, then computes their average.
