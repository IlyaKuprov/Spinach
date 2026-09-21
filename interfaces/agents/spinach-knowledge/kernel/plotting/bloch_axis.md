# kernel/plotting/bloch_axis.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/bloch_axis.m`
- Signature: `[ax,ay,az]=bloch_axis(x,y,z)`
- Total lines: 64

## Purpose

Reconstructs the instantaneous Bloch equation rotation axis of from a 3D magnetisation trajectory. Syntax: [ax,ay,az]=bloch_axis(x,y,z)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `fdvec()`, `ismatrix()`, `any()`, `isequal()`.
