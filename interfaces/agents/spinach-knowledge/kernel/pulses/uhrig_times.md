# kernel/pulses/uhrig_times.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/pulses/uhrig_times.m`
- Signature: `time_delays=uhrig_times(T,N)`
- Total lines: 61

## Purpose

Uhrig's UDD decoupling sequence timings. Syntax: time_delays=uhrig_times(T,N)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 25-26: Check consistency; implemented by `grumble(T,N)`.
- Lines 28-29: Use the formula from WSW's 2009 JCP paper; implemented by `time_positions=T*(sin(pi*(1:N)/(2*N+2)).^2-0.5)`.
- Lines 31-32: Convert positions to delays; implemented by `time_delays=diff(time_positions)`.
- Lines 34-35: Add the starting and the trailing delay; implemented by `chunk=(T-sum(time_delays))/2`.

### Key state/data transformations

- Lines 29: computes `time_positions` using `time_positions=T*(sin(pi*(1:N)/(2*N+2)).^2-0.5)`.
- Lines 32: computes `time_delays` using `time_delays=diff(time_positions)`.
- Lines 35: computes `chunk` using `chunk=(T-sum(time_delays))/2`.

### Local helper functions

- Line 41: `grumble()` — `function grumble(T,N)`. "In a legislative effort to force vendors to hand over the plaintext contents of encrypted communications, Prime Minister Malcom Turnbull
  - Representative operation: `if (~isnumeric(T))||(~isreal(T))||(~isscalar(T))||(~isfinite(T))||(T<=0)`.
  - Representative operation: `error('T must be a finite positive real scalar.')`.

## Parameters / inputs

- T -total duration of the sequence (sum of all delays)
- N -number of pulses in the sequence

## Outputs

- time_delays -list of delays between ideal pulses in
- the UDD sequence; the first pulse goes
- after the first delay, and there is a
- delay after the last pulse

## Implementation structure

- Uhrig's UDD decoupling sequence timings. Syntax:
- time_delays=uhrig_times(T,N)
- T -total duration of the sequence (sum of all delays)
- N -number of pulses in the sequence
- time_delays -list of delays between ideal pulses in
- the UDD sequence; the first pulse goes
- after the first delay, and there is a
- delay after the last pulse
- Check consistency
- Use the formula from WSW's 2009 JCP paper
- Convert positions to delays
- Add the starting and the trailing delay

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `diff()`, `isscalar()`.
