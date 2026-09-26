# kernel/conventions/transforms/ang2cgsppm.m

- Signature: `cgsppm=ang2cgsppm(ang)`

## Purpose

Converts magnetic susceptibility from the Angstrom^3 units required by Spinach pseudocontact shift functionality into the cgs-ppm (aka cm^3/mol) units quoted by quantum chemist- ry packages. Syntax: cgsppm=ang2cgsppm(ang)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

## Parameters / inputs

- ang -an array of values in cubic Angstrom

## Outputs

- cgsppm -an array of values in cgs-ppm

## Implementation structure

- Converts magnetic susceptibility from the Angstrom^3 units
- required by Spinach pseudocontact shift functionality into
- the cgs-ppm (aka cm^3/mol) units quoted by quantum chemist-
- ry packages. Syntax:
- cgsppm=ang2cgsppm(ang)
- ang -an array of values in cubic Angstrom
- cgsppm -an array of values in cgs-ppm
- Check consistency
- Do the calculation
- Consistency enforcement
- No artist tolerates reality.
- Friedrich Nietzsche
