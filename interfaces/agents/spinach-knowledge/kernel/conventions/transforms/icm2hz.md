# kernel/conventions/transforms/icm2hz.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/icm2hz.m`
- Signature: `hz=icm2hz(icm)`
- Total lines: 39

## Purpose

Converts cm^-1 units used in spectroscopy into Hz units preferred in magnetic resonance. Syntax: hz=icm2hz(icm) Arrays of any dimensions are supported. Parameters: icm -an array of values in inverse centimetres

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 20-21: Check consistency; implemented by `grumble(icm)`.
- Lines 23-24: Run the conversion; implemented by `hz=100*299792458*icm`.

### Key state/data transformations

- Lines 24: computes `hz` using `hz=100*299792458*icm`.

### Local helper functions

- Line 29: `grumble()` — `function grumble(icm)`. Satire thrives where the usual checks on human folly fail.
  - Representative operation: `if (~isnumeric(icm))||(~isreal(icm))`.
  - Representative operation: `error('the argument must be an array of real numbers.')`.

## Outputs

- hz -an array of values in Hz

## Implementation structure

- Converts cm^-1 units used in spectroscopy into Hz units
- preferred in magnetic resonance. Syntax:
- hz=icm2hz(icm)
- Arrays of any dimensions are supported. Parameters:
- icm -an array of values in inverse centimetres
- hz -an array of values in Hz
- Check consistency
- Run the conversion
- Consistency enforcement
- Satire thrives where the usual checks on human
- folly fail.
- The Economist

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`.
