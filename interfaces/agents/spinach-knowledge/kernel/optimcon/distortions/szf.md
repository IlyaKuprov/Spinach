# kernel/optimcon/distortions/szf.m

- Signature: `[w,J]=szf(w,z)`

## Purpose

Applies a discrete single-zero filter: Y(k)=X(k)/(1-z)-z*X(k-1)/(1-z); to a Spinach optimal control module waveform. Treats odd rows of multi-row waveform arrays as real, and even rows as imaginary, components of a complex signal. Syntax: [w,J]=szf(w,z)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

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
