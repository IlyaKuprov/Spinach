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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 30-31: Check consistency; implemented by `grumble(npoints,extents)`.
- Lines 33-34: Decide the dimensionality; implemented by `switch numel(npoints)`.
- Lines 38-39: Get differentiation matrices; implemented by `[~,Dxx]=fourdif(npoints(1),2)`.
- Lines 41-42: Normalize differentiation matrices; implemented by `Dxx=(2*pi/extents(1))^2*Dxx`.
- Lines 44-45: Compute the Laplacian; implemented by `L=Dxx`.
- Lines 57-59: Compute the Laplacian; implemented by `L=kron(Dyy,speye(npoints(1)))+ kron(speye(npoints(2)),Dxx)`.
- Lines 73-76: Compute the Laplacian; implemented by `L=kron(kron(Dzz,speye(npoints(2))),speye(npoints(1)))+ kron(kron(speye(npoints(3)),Dyy),speye(npoints(1)))+ kron(kron(speye(npoints(3)),speye(npoints(2))),Dxx)`.
- Lines 80-81: Complain and bomb out; implemented by `error('incorrect number of spatial dimensions.')`.

### Control flow inferred from the code

- Line 34: dispatches on `numel(npoints)`; cases `1`, `2`, `3`.

### Key state/data transformations

- Lines 39: computes `[~,Dxx]` using `[~,Dxx]=fourdif(npoints(1),2)`.
- Lines 42: computes `Dxx` using `Dxx=(2*pi/extents(1))^2*Dxx`.
- Lines 45: computes `L` using `L=Dxx`.
- Lines 51: computes `[~,Dyy]` using `[~,Dyy]=fourdif(npoints(2),2)`.
- Lines 55: computes `Dyy` using `Dyy=(2*pi/extents(2))^2*Dyy`.
- Lines 66: computes `[~,Dzz]` using `[~,Dzz]=fourdif(npoints(3),2)`.
- Lines 71: computes `Dzz` using `Dzz=(2*pi/extents(3))^2*Dzz`.

### Local helper functions

- Line 88: `grumble()` — `function grumble(npoints,extents)`. Of course it's the same old story. Truth usually is the same old story.
  - Representative operation: `if (~isnumeric(npoints))||(~isreal(npoints))||(any(npoints<1))||any(mod(npoints,1)~=0)`.
  - Representative operation: `error('npoints must be a three-element vector of positive integers.')`.

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
