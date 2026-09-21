# kernel/optimcon/distortions/spf.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/distortions/spf.m`
- Signature: `[w,J]=spf(w,p)`
- Total lines: 125

## Purpose

Applies a discrete single-pole filter: Y(n)=(1-p)*X(n)+p*Y(n-1) to a Spinach optimal control module waveform. Treats odd rows of multi-row waveform arrays as real, and even rows as imaginary, components of a complex signal. Syntax: [w,J]=spf(w,p)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `distort()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- w -waveform, one time slice per column, and
- rows arranged as XYXY... with respect to
- in-phase and quadrature parts on each
- control channel
- p -a vector (one element per XY control pair)
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

- Applies a discrete single-pole filter:
- Y(n)=(1-p)*X(n)+p*Y(n-1)
- to a Spinach optimal control module waveform. Treats odd
- rows of multi-row waveform arrays as real, and even rows
- as imaginary, components of a complex signal. Syntax:
- [w,J]=spf(w,p)
- w -waveform, one time slice per column, and
- rows arranged as XYXY... with respect to
- in-phase and quadrature parts on each
- control channel
- p -a vector (one element per XY control pair)
- containing the filter coefficient:

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `distort()`, `dlfeval()`, `dlarray()`, `extractdata()`, `dims()`, `inp()`, `transpose()`, `w_dist()`, `dljacobian()`, `isvector()`, `any()`.
