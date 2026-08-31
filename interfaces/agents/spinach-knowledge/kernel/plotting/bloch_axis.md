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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 22-23: Check consistency; implemented by `grumble(x,y,z)`.
- Lines 25-26: Get first and second derivatives; implemented by `dx_dt=fdvec(x,5,1); d2x_dt2=fdvec(x,5,2)`.
- Lines 30-31: Get instantaneous rotation axis; implemented by `ax=dy_dt.*d2z_dt2-dz_dt.*d2y_dt2`.

### Key state/data transformations

- Lines 26: computes `dx_dt` using `dx_dt=fdvec(x,5,1); d2x_dt2=fdvec(x,5,2)`.
- Lines 27: computes `dy_dt` using `dy_dt=fdvec(y,5,1); d2y_dt2=fdvec(y,5,2)`.
- Lines 28: computes `dz_dt` using `dz_dt=fdvec(z,5,1); d2z_dt2=fdvec(z,5,2)`.
- Lines 31: computes `ax` using `ax=dy_dt.*d2z_dt2-dz_dt.*d2y_dt2`.
- Lines 32: computes `ay` using `ay=dz_dt.*d2x_dt2-dx_dt.*d2z_dt2`.
- Lines 33: computes `az` using `az=dx_dt.*d2y_dt2-dy_dt.*d2x_dt2`.

### Local helper functions

- Line 38: `grumble()` — `function grumble(x,y,z)`.
  - Representative operation: `if (~isnumeric(x))||(~isreal(x))|| (~ismatrix(x))||any(~isfinite(x),'all')`.
  - Representative operation: `(~ismatrix(x))||any(~isfinite(x),'all')`.

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
