# kernel/conventions/transforms/mev2hz.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/mev2hz.m`
- Signature: `hz=mev2hz(mev)`
- Total lines: 42

## Purpose

Converts meV energy units used in solid state physics and phonon spectroscopy into Hz units preferred in magnetic resonance. Syntax: hz=mev2hz(mev) Arrays of any dimensions are supported. Parameters: mev -an array of values in milli-electronvolts

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 21-22: Check consistency; implemented by `grumble(mev)`.
- Lines 24-25: Run the conversion; implemented by `hz=1e-3*1.602176634e-19*mev/6.62607015e-34`.

### Key state/data transformations

- Lines 25: computes `hz` using `hz=1e-3*1.602176634e-19*mev/6.62607015e-34`.

### Local helper functions

- Line 30: `grumble()` — `function grumble(mev)`. O God, I could be bounded in a nutshell, and count myself a king of infinite space, were it not that I
  - Representative operation: `if (~isnumeric(mev))||(~isreal(mev))`.
  - Representative operation: `error('the argument must be an array of real numbers.')`.

## Outputs

- hz -an array of values in Hz

## Implementation structure

- Converts meV energy units used in solid state physics and
- phonon spectroscopy into Hz units preferred in magnetic
- resonance. Syntax:
- hz=mev2hz(mev)
- Arrays of any dimensions are supported. Parameters:
- mev -an array of values in milli-electronvolts
- hz -an array of values in Hz
- Check consistency
- Run the conversion
- Consistency enforcement
- O God, I could be bounded in a nutshell, and count
- myself a king of infinite space, were it not that I

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`.
