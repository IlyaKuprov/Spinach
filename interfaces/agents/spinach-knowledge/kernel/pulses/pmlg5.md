# kernel/pulses/pmlg5.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/pulses/pmlg5.m`
- Signature: `phi=pmlg5(n)`
- Total lines: 48

## Purpose

PMLG5 phase sequence as described in the paper by Vinogradova, Madhu and Vega (https://doi.org/10.1016/S0009-2614(99)01174-4).

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 22-23: Check consistency; implemented by `grumble(n)`.
- Lines 25-29: PMLG5 phase sequence; implemented by `phi_sequence=[339.22 297.65 256.08 214.51 172.94 352.94 34.51 76.08 117.65 159.22 159.22 117.65 76.08 34.51 352.94 172.94 214.51 256.08 297.65 339.22]`.
- Lines 31-32: Loop correctly over; implemented by `phi=(pi/180)*phi_sequence(mod(n-1,20)+1)`.

### Key state/data transformations

- Lines 26-29: computes `phi_sequence` using `phi_sequence=[339.22 297.65 256.08 214.51 172.94 352.94 34.51 76.08 117.65 159.22 159.22 117.65 76.08 34.51 352.94 172.94 214.51 256.08 297.65 339.22]`.
- Lines 32: computes `phi` using `phi=(pi/180)*phi_sequence(mod(n-1,20)+1)`.

### Local helper functions

- Line 37: `grumble()` — `function grumble(n)`. One man's crappy software is another man's full time job.
  - Representative operation: `if (~isnumeric(n))||(~isreal(n))|| (~isscalar(n))||(n<1)||(mod(n,1)~=0)`.
  - Representative operation: `(~isscalar(n))||(n<1)||(mod(n,1)~=0)`.

## Syntax

```matlab
phi=spinal(n)
```

## Parameters / inputs

- n -a positive integer number

## Outputs

- phi -the phase of the n-th pulse in
- PMLG sequence, radians

## Implementation structure

- PMLG5 phase sequence as described in the paper by Vinogradova,
- Madhu and Vega (https://doi.org/10.1016/S0009-2614(99)01174-4).
- phi=spinal(n)
- n -a positive integer number
- phi -the phase of the n-th pulse in
- PMLG sequence, radians
- Check consistency
- PMLG5 phase sequence
- Loop correctly over
- Consistency enforcement
- One man's crappy software is another
- man's full time job.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `phi_sequence()`, `isscalar()`.
