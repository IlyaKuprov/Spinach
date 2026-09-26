# kernel/optimcon/distortions/firf.m

- Signature: `[w,J]=firf(w,ker)`

## Purpose

Applies an FIR convolution filter to a Spinach optimal control module waveform. Treats odd rows of multi-row waveform arrays as real, and even rows as imaginary, components of a complex signal. The distal end of the convolution is truncated so the output has the same number of samples as the input. Syntax: [w,J]=firf(w,ker)

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
- ker -a vector of FIR filter coefficients

## Outputs

- w -distorted waveform, same dimension as the
- input waveform; leaving sufficient ring-
- down margin is the user's responsibility
- J -Jacobian matrix with respect to vectorisa-
- tions of the output and the input arrays

## Implementation structure

- Applies an FIR convolution filter to a Spinach optimal control
- module waveform. Treats odd rows of multi-row waveform arrays
- as real, and even rows as imaginary, components of a complex
- signal. The distal end of the convolution is truncated so the
- output has the same number of samples as the input. Syntax:
- [w,J]=firf(w,ker)
- w -waveform, one time slice per column, and
- rows arranged as XYXY... with respect to
- in-phase and quadrature parts on each
- control channel
- ker -a vector of FIR filter coefficients
- w -distorted waveform, same dimension as the
