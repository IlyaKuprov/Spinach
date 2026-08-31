# kernel/conventions/transforms/hz2ppm.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/hz2ppm.m`
- Signature: `ppm=hz2ppm(hz,B0,nucleus)`
- Total lines: 55

## Purpose

Converts resonance offsets into chemical shifts. Syntax: ppm=hz2ppm(hz,B0,nucleus)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 26-27: Check consistency; implemented by `grumble(hz,B0,nucleus)`.
- Lines 29-30: Calculate chemical shift in Hz; implemented by `ppm=1e6*(2*pi*hz)/(B0*spin(nucleus))`.

### Key state/data transformations

- Lines 30: computes `ppm` using `ppm=1e6*(2*pi*hz)/(B0*spin(nucleus))`.

### Local helper functions

- Line 35: `grumble()` — `function grumble(hz,B0,nucleus)`.
  - Representative operation: `if (~isnumeric(hz))||(~isreal(hz))`.
  - Representative operation: `error('resonance offset must be real.')`.

## Parameters / inputs

- hz -resonance offset in Hz
- B0 -magnet induction, Tesla
- nucleus -a string specifying the isotope, e.g. '1H'

## Outputs

- ppm -chemical shift in ppm
- Note: signs of the magnetogyric ratios are preserved.

## Implementation structure

- Converts resonance offsets into chemical shifts. Syntax:
- ppm=hz2ppm(hz,B0,nucleus)
- hz -resonance offset in Hz
- B0 -magnet induction, Tesla
- nucleus -a string specifying the isotope, e.g. '1H'
- ppm -chemical shift in ppm
- Note: signs of the magnetogyric ratios are preserved.
- Check consistency
- Calculate chemical shift in Hz
- Consistency enforcement
- Somebody Else's Wife: oh, I am so tired of being boringly taken
- for granted, he never really acknowledges

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `ischar()`.
