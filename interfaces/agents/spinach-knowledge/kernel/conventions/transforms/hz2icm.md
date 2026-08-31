# kernel/conventions/transforms/hz2icm.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/hz2icm.m`
- Signature: `icm=hz2icm(hz)`
- Total lines: 40

## Purpose

Converts Hz units used in magnetic resonance into cm^-1 units used in spectroscopy. Syntax: icm=hz2icm(hz) Arrays of any dimensions are supported. Parameters: hz -an array of values in Hz

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 20-21: Check consistency; implemented by `grumble(hz)`.
- Lines 23-24: Run the conversion; implemented by `icm=hz/(100*299792458)`.

### Key state/data transformations

- Lines 24: computes `icm` using `icm=hz/(100*299792458)`.

### Local helper functions

- Line 29: `grumble()` — `function grumble(hz)`. "Damn, would I have to be nice to everyone for two years?!" IK, upon being informed that he was to
  - Representative operation: `if (~isnumeric(hz))||(~isreal(hz))`.
  - Representative operation: `error('the argument must be an array of real numbers.')`.

## Outputs

- icm -an array of values in inverse centimetres

## Implementation structure

- Converts Hz units used in magnetic resonance into cm^-1 units
- used in spectroscopy. Syntax:
- icm=hz2icm(hz)
- Arrays of any dimensions are supported. Parameters:
- hz -an array of values in Hz
- icm -an array of values in inverse centimetres
- Check consistency
- Run the conversion
- Consistency enforcement
- "Damn, would I have to be nice to everyone for two years?!"
- IK, upon being informed that he was to
- organise the 48th ESR Group Conference

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`.
