# kernel/optimcon/distortions/szf.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/distortions/szf.m`
- Signature: `[w,J]=szf(w,z)`
- Total lines: 131

## Purpose

Applies a discrete single-zero filter: Y(k)=X(k)/(1-z)-z*X(k-1)/(1-z); to a Spinach optimal control module waveform. Treats odd rows of multi-row waveform arrays as real, and even rows as imaginary, components of a complex signal. Syntax: [w,J]=szf(w,z)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `distort()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 44-45: Check consistency; implemented by `grumble(w,z)`.
- Lines 47-48: Autodiff wrapper; implemented by `if nargout<2`.
- Lines 50-51: Plain call; implemented by `w=distort(w(:),z,size(w))`.
- Lines 55-56: Autodiff call including Jacobian; implemented by `[w,J]=dlfeval(@distort,dlarray(w(:)),z,size(w))`.
- Lines 58-59: Strip autodiff rigging; kill Wirtinger terms; implemented by `w=extractdata(w); J=extractdata(J); J=real(J)`.

### Control flow inferred from the code

- Line 48: conditional branch on `nargout<2`.

### Key state/data transformations

- Lines 51: computes `w` using `w=distort(w(:),z,size(w))`.
- Lines 56: computes `[w,J]` using `[w,J]=dlfeval(@distort,dlarray(w(:)),z,size(w))`.

### Local helper functions

- Line 66: `distort()` — `function [w_dist,J]=distort(w,z,dims)`. Fold into physical dimensions
  - Representative operation: `inp=reshape(w,dims); nrows=dims(1)`.
  - Representative operation: `ncols=dims(2); nchannels=nrows/2`.
- Line 102: `grumble()` — `function grumble(w,z)`.
  - Representative operation: `if (~isnumeric(w))||(~isreal(w))`.
  - Representative operation: `error('w must be an array of real numbers.')`.

## Parameters / inputs

- w -waveform, one time slice per column, and
- rows arranged as XYXY... with respect to
- in-phase and quadrature parts on each
- control channel
- z -a vector (one element per XY control pair)
- containing the filter coefficient:
- p=exp(-r*dt+1i*(omega-omega_rf)*dt)
- where r is th damping rate, omega is the
- pole frequency, omega_rf is the rotating
- frame frequency, and dt is the time dis-
- cretisation step.

## Outputs

- w -distorted waveform, same dimension as the
- input waveform; leaving sufficient ring-
- down margin is the user's responsibility
- J -Jacobian matrix with respect to vectorisa-
- tions of the output and the input arrays

## Implementation structure

- Applies a discrete single-zero filter:
- Y(k)=X(k)/(1-z)-z*X(k-1)/(1-z);
- to a Spinach optimal control module waveform. Treats odd
- rows of multi-row waveform arrays as real, and even rows
- as imaginary, components of a complex signal. Syntax:
- [w,J]=szf(w,z)
- w -waveform, one time slice per column, and
- rows arranged as XYXY... with respect to
- in-phase and quadrature parts on each
- control channel
- z -a vector (one element per XY control pair)
- containing the filter coefficient:

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `distort()`, `dlfeval()`, `dlarray()`, `extractdata()`, `dims()`, `inp()`, `transpose()`, `w_dist()`, `dljacobian()`, `isvector()`, `any()`.
