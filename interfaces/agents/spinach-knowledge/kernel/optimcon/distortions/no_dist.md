# kernel/optimcon/distortions/no_dist.m

- Signature: `[w,J]=no_dist(w)`

## Purpose

A distortion function that applies no distortion and therefore has a unit Jacobian. Syntax: [w,J]=no_dist(w)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.

## Numerical / algorithmic content

## Parameters / inputs

- w -waveform, a numerical array

## Outputs

- w -the same waveform as the input
- J -a sparse unit matrix with the dimension mat-
- ching the vectorisation of the input

## Implementation structure

- A distortion function that applies no distortion and therefore
- has a unit Jacobian. Syntax:
- [w,J]=no_dist(w)
- w -waveform, a numerical array
- w -the same waveform as the input
- J -a sparse unit matrix with the dimension mat-
- ching the vectorisation of the input
- Check consistency
- Return a unit Jacobian if asked
- Consistency enforcement
- If I only knew how I could get mathematicians interested in
- transformation groups and the treatment of differential equ-
