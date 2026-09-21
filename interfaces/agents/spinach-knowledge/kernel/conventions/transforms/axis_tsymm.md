# kernel/conventions/transforms/axis_tsymm.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/axis_tsymm.m`
- Signature: `A=axis_tsymm(T,a)`
- Total lines: 61

## Purpose

Roughly averages an interaction tensor with respect to the rotation around a user-specified axis. Syntax: T=axis_tsymm(T,a)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `anax2dcm()`.
