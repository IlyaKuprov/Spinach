# kernel/conventions/transforms/ang2cgsppm.m

- Signature: `cgsppm=ang2cgsppm(ang)`

## Purpose

Converts magnetic susceptibility from the cubic-angstrom units required by Spinach pseudocontact-shift functionality to the cgs-ppm (cm^3/mol) units quoted by quantum chemistry packages.

## Physical / mathematical content

- Conversion: `cgsppm = 6.02214129e23 * ang / (4*pi*1e18)`.

## Numerical / algorithmic content

- The calculation applies the conversion factor to the input array.
- The input must be numeric; otherwise, the function raises an error.

## Parameters / inputs

- ang -an array of values in cubic Angstrom

## Outputs

- cgsppm -an array of values in cgs-ppm

## Implementation structure

- Checks that `ang` is numeric, then calculates `cgsppm` using the conversion formula.
