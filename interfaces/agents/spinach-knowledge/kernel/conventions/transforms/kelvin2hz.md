# kernel/conventions/transforms/kelvin2hz.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/kelvin2hz.m`
- Signature: `hz=kelvin2hz(kelvin)`
- Total lines: 44

## Purpose

Converts Kelvin energy units used for Debye temperatures and thermal energy scales in solid state physics into Hz units preferred in magnetic resonance. Syntax: hz=kelvin2hz(kelvin) Arrays of any dimensions are supported. Parameters: kelvin -an array of values in Kelvin

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 21-22: Check consistency; implemented by `grumble(kelvin)`.
- Lines 24-25: Run the conversion; implemented by `hz=1.380649e-23*kelvin/6.62607015e-34`.

### Key state/data transformations

- Lines 25: computes `hz` using `hz=1.380649e-23*kelvin/6.62607015e-34`.

### Local helper functions

- Line 30: `grumble()` — `function grumble(kelvin)`. When you can measure what you are speaking about, and express it in numbers, you know something about it; but
  - Representative operation: `if (~isnumeric(kelvin))||(~isreal(kelvin))`.
  - Representative operation: `error('the argument must be an array of real numbers.')`.

## Outputs

- hz -an array of values in Hz

## Implementation structure

- Converts Kelvin energy units used for Debye temperatures and
- thermal energy scales in solid state physics into Hz units
- preferred in magnetic resonance. Syntax:
- hz=kelvin2hz(kelvin)
- Arrays of any dimensions are supported. Parameters:
- kelvin -an array of values in Kelvin
- hz -an array of values in Hz
- Check consistency
- Run the conversion
- Consistency enforcement
- When you can measure what you are speaking about, and
- express it in numbers, you know something about it; but

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`.
