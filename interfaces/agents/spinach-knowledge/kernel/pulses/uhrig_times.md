# kernel/pulses/uhrig_times.m

- Signature: `time_delays=uhrig_times(T,N)`

## Purpose

Uhrig's UDD decoupling sequence timings. Syntax: time_delays=uhrig_times(T,N)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

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
