# kernel/optimcon/distortions/non_orth.m

- Signature: `[w,J]=non_orth(w,xy_ang)`

## Purpose

Non-orthogonal channel distortion model. Treats odd rows of multi-row waveform arrays as in-phase channels, and even rows as out-of-phase channels. The in-phase channel is kept fixed; the out-of-phase channel is tilted so that its true angle to the in-phase channel is user-specified. Syntax: [w,J]=non_orth(w,xy_ang)

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
- xy_ang -true angle, in degrees, between the instru-
- ment output directions of each X,Y control
- pair; may be a scalar or one value per pair,
- with 90 degrees corresponding to no distortion

## Outputs

- w -distorted waveform, same dimension as the
- input waveform
- J -Jacobian matrix with respect to vectorisa-
- tions of the output and the input arrays

## Implementation structure

- Non-orthogonal channel distortion model. Treats odd rows of
- multi-row waveform arrays as in-phase channels, and even rows
- as out-of-phase channels. The in-phase channel is kept fixed;
- the out-of-phase channel is tilted so that its true angle to
- the in-phase channel is user-specified. Syntax:
- [w,J]=non_orth(w,xy_ang)
- w -waveform, one time slice per column, and
- rows arranged as XYXY... with respect to
- in-phase and quadrature parts on each
- control channel
- xy_ang -true angle, in degrees, between the instru-
- ment output directions of each X,Y control
