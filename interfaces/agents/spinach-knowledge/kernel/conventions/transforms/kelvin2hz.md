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
