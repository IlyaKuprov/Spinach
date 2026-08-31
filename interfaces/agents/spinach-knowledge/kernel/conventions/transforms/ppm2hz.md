# kernel/conventions/transforms/ppm2hz.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/ppm2hz.m`
- Signature: `hz=ppm2hz(ppm,B0,nucleus)`
- Total lines: 50

## Purpose

Converts chemical shifts into resonance offsets. Syntax: hz=ppm2hz(ppm,B0,nucleus)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 26-27: Check consistency; implemented by `grumble(ppm,B0,nucleus)`.
- Lines 29-30: Calculate chemical shift in Hz; implemented by `hz=1e-6*ppm*(B0*spin(nucleus)/(2*pi))`.

### Key state/data transformations

- Lines 30: computes `hz` using `hz=1e-6*ppm*(B0*spin(nucleus)/(2*pi))`.

### Local helper functions

- Line 35: `grumble()` — `function grumble(ppm,B0,nucleus)`.
  - Representative operation: `if (~isnumeric(ppm))||(~isreal(ppm))`.
  - Representative operation: `error('chemical shift must be real.')`.

## Parameters / inputs

- ppm -chemical shift in ppm
- B0 -magnet induction, Tesla
- nucleus -a string specifying the isotope, e.g. '1H'

## Outputs

- hz -resonance offset in Hz
- Note: signs of the magnetogyric ratios are preserved.

## Implementation structure

- Converts chemical shifts into resonance offsets. Syntax:
- hz=ppm2hz(ppm,B0,nucleus)
- ppm -chemical shift in ppm
- B0 -magnet induction, Tesla
- nucleus -a string specifying the isotope, e.g. '1H'
- hz -resonance offset in Hz
- Note: signs of the magnetogyric ratios are preserved.
- Check consistency
- Calculate chemical shift in Hz
- Consistency enforcement
- "Smoking -NO HYDROGEN!"
- Safety warning on Anatole Abragam's door

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `ischar()`.
