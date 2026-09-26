# kernel/plotting/bloch_axis.m

- Signature: `[ax,ay,az]=bloch_axis(x,y,z)`

## Purpose

Reconstructs the instantaneous Bloch equation rotation axis of from a 3D magnetisation trajectory. Syntax: [ax,ay,az]=bloch_axis(x,y,z)

## Physical / mathematical content

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.

## Parameters / inputs

- x, y, z -row vectors of equal length con-
- taining the trajectory
- Output:
- ax, ay, az -row vectors of equal length con-
- taining the instantaneous axis

## Implementation structure

- Reconstructs the instantaneous Bloch equation rotation axis of
- from a 3D magnetisation trajectory. Syntax:
- [ax,ay,az]=bloch_axis(x,y,z)
- x, y, z -row vectors of equal length con-
- taining the trajectory
- Output:
- ax, ay, az -row vectors of equal length con-
- taining the instantaneous axis
- Check consistency
- Get first and second derivatives
- Get instantaneous rotation axis
- Consistency enforcement
