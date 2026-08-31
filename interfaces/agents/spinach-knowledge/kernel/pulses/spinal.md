# kernel/pulses/spinal.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/pulses/spinal.m`
- Signature: `phi=spinal(n)`
- Total lines: 50

## Purpose

SPINAL phase sequences as described in the paper by Fung, Khitrin and Ermolaev (https://doi.org/10.1006/jmre.1999.1896). Syntax: phi=spinal(n)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 22-23: Check consistency; implemented by `grumble(n)`.
- Lines 25-31: Spinal phase sequence; implemented by `phi_sequence=[10 350 15 345 20 340 15 345 350 10 345 15 340 20 345 15 350 10 345 15 340 20 345 15 10 350 15 345 20 340 15 345 350 10 345 15 340 20 345 15 10 350 15 345 2…`.
- Lines 33-34: Loop correctly over; implemented by `phi=(pi/180)*phi_sequence(mod(n-1,64)+1)`.

### Key state/data transformations

- Lines 26-31: computes `phi_sequence` using `phi_sequence=[10 350 15 345 20 340 15 345 350 10 345 15 340 20 345 15 350 10 345 15 340 20 345 15 10 350 15 345 20 340 15 345 350 10 345 15 340 20 345 15 10 350 15 345 2…`.
- Lines 34: computes `phi` using `phi=(pi/180)*phi_sequence(mod(n-1,64)+1)`.

### Local helper functions

- Line 39: `grumble()` — `function grumble(n)`. The key to performance is elegance, not battalions of special cases.
  - Representative operation: `if (~isnumeric(n))||(~isreal(n))|| (~isscalar(n))||(n<1)||(mod(n,1)~=0)`.
  - Representative operation: `(~isscalar(n))||(n<1)||(mod(n,1)~=0)`.

## Parameters / inputs

- n -a positive integer number

## Outputs

- phi -the phase of the n-th pulse in
- SPINAL sequence, radians

## Implementation structure

- SPINAL phase sequences as described in the paper by Fung, Khitrin
- and Ermolaev (https://doi.org/10.1006/jmre.1999.1896). Syntax:
- phi=spinal(n)
- n -a positive integer number
- phi -the phase of the n-th pulse in
- SPINAL sequence, radians
- Check consistency
- Spinal phase sequence
- Loop correctly over
- Consistency enforcement
- The key to performance is elegance, not
- battalions of special cases.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `phi_sequence()`, `isscalar()`.
