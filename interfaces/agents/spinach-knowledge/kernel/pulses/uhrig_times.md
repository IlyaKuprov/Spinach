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
