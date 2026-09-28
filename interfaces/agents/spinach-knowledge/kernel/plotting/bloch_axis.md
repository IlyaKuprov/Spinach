# kernel/plotting/bloch_axis.m

- Signature: `[ax,ay,az]=bloch_axis(x,y,z)`

## Purpose

Reconstructs the instantaneous Bloch equation rotation axis of from a 3D magnetisation trajectory. Syntax: [ax,ay,az]=bloch_axis(x,y,z)

## Physical / mathematical content

For a 3D magnetisation trajectory, returns the instantaneous rotation-axis vector from the cross product of the first and second derivatives.

## Numerical / algorithmic content

Each component’s first and second derivatives are computed with `fdvec(component,5,1)` and `fdvec(component,5,2)`. The returned components are `ax=dy.*d2z-dz.*d2y`, `ay=dz.*d2x-dx.*d2z`, and `az=dx.*d2y-dy.*d2x`.

## Parameters / inputs

- `x`, `y`, `z` — equal-length row vectors containing the trajectory; the source checks that they are real, finite, numeric, and have identical dimensions.

## Outputs

- `ax`, `ay`, `az` — components of the instantaneous rotation axis, with dimensions matching the input vectors.
