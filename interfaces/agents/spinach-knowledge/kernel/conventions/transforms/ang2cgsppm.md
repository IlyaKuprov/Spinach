# kernel/conventions/transforms/ang2cgsppm.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/ang2cgsppm.m`
- Signature: `cgsppm=ang2cgsppm(ang)`
- Total lines: 40

## Purpose

Converts magnetic susceptibility from the Angstrom^3 units required by Spinach pseudocontact shift functionality into the cgs-ppm (aka cm^3/mol) units quoted by quantum chemist- ry packages. Syntax: cgsppm=ang2cgsppm(ang)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 22-23: Check consistency; implemented by `grumble(ang)`.
- Lines 25-26: Do the calculation; implemented by `cgsppm=6.02214129e23*ang/(4*pi*1e18)`.

### Key state/data transformations

- Lines 26: computes `cgsppm` using `cgsppm=6.02214129e23*ang/(4*pi*1e18)`.

### Local helper functions

- Line 31: `grumble()` — `function grumble(ang)`. No artist tolerates reality. Friedrich Nietzsche
  - Representative operation: `if ~isnumeric(ang)`.
  - Representative operation: `error('input must be numeric.')`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`.
