# kernel/conventions/transforms/g2freq.m

- Signature: `f=g2freq(g,B)`

## Purpose

Converts g-tensor units into electron Zeeman frequency units. Syntax: f=g2freq(g,B)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

## Parameters / inputs

- g -g-values, scalar or array
- B -magnetic field in Tesla

## Outputs

- f -frequency in Hz

## Implementation structure

- Converts g-tensor units into electron Zeeman frequency
- units. Syntax:
- f=g2freq(g,B)
- g - g-values, scalar or array
- B - magnetic field in Tesla
- f - frequency in Hz
- Check consistency
- Get the free electron carrier frequency
- Scale the frequency
- Consistency enforcement
- "Authors are listed in order of degree of belief in
- the central thesis."
