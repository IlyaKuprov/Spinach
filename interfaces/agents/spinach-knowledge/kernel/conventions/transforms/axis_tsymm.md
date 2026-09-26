# kernel/conventions/transforms/axis_tsymm.m

- Signature: `A=axis_tsymm(T,a)`

## Purpose

Roughly averages an interaction tensor with respect to the rotation around a user-specified axis. Syntax: T=axis_tsymm(T,a)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

## Parameters / inputs

- T -3x3 real interaction tensor
- a -3x1 real vector specifying the rotation axis

## Outputs

- A -3x3 real interaction tensor averaged over
- the rotation around the specified axis

## Implementation structure

- Roughly averages an interaction tensor with respect to the
- rotation around a user-specified axis. Syntax:
- T=axis_tsymm(T,a)
- T -3x3 real interaction tensor
- a -3x1 real vector specifying the rotation axis
- A -3x3 real interaction tensor averaged over
- the rotation around the specified axis
- Check consistency
- Preallocate average
- Loop over the full rotation
- Compute the average
- Consistency enforcement
