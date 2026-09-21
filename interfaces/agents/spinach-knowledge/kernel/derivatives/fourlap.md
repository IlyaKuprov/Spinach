# kernel/derivatives/fourlap.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/derivatives/fourlap.m`
- Signature: `L=fourlap(npoints,extents)`
- Total lines: 101

## Purpose

Returns a Fourier spectral representation of the Laplacian acting on a 3D data array. Syntax: L=fourlap(npoints,extents)

## Physical / mathematical content

- Derivative utilities. These routines compute finite-difference, analytical, or optimisation-oriented derivatives needed for sensitivity analysis, fitting, and optimal control.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- npoints -a three-element vector specifying the number of
- discretization points in each dimension of the
- 3D cube of data that the operator will be acting
- on, ordered as [X Y Z].
- extents -a three-element vector specifying axis extents,
- ordered as [X Y Z].

## Outputs

- L -Fourier spectral Laplacian, a sparse matrix designed to act
- on the vectorization of the 3D data array. The dimensions of
- the data array are assumed to be ordered as [X Y Z].
- Note: periodic boundary conditions.

## Implementation structure

- Returns a Fourier spectral representation of the Laplacian acting
- on a 3D data array. Syntax:
- L=fourlap(npoints,extents)
- npoints - a three-element vector specifying the number of
- discretization points in each dimension of the
- 3D cube of data that the operator will be acting
- on, ordered as [X Y Z].
- extents - a three-element vector specifying axis extents,
- ordered as [X Y Z].
- L -Fourier spectral Laplacian, a sparse matrix designed to act
- on the vectorization of the 3D data array. The dimensions of
- the data array are assumed to be ordered as [X Y Z].

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `fourdif()`, `npoints()`, `extents()`, `speye()`, `any()`.
