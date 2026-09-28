# kernel/conventions/transforms/cgsppm2ang.m

- Signature: `ang=cgsppm2ang(cgsppm)`

## Purpose

Converts magnetic susceptibility from the cgs-ppm (aka cm^3/mol) units quoted by quantum chemistry packages into Angstrom^3 units required by Spinach pseudocontact shift functionality. Syntax: ang=cgsppm2ang(cgsppm)

## Physical / mathematical content
Converts susceptibility values from cgs-ppm (described in the source as cm^3/mol) to cubic angstroms using `ang = (4*pi*1e18/6.02214129e23)*cgsppm`. The scalar factor is applied elementwise, preserving the input array's size.

## Numerical / algorithmic content

## Parameters / inputs

- cgsppm -any numerical array of susceptibility
- values in cgs-ppm

## Outputs

- ang -array of the same size with suscepti-
- bility values in cubic Angstrom

## Implementation structure
The function calls the local `grumble` helper, which raises `input must be numeric.` if `cgsppm` is not numeric. It then computes `ang=4*pi*1e18*cgsppm/6.02214129e23;`. No other input checks or transformations are performed.
