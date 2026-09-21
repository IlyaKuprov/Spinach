# kernel/pulses/pulse_shape.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/pulses/pulse_shape.m`
- Signature: `waveform=pulse_shape(pulse_name,npoints)`
- Total lines: 74

## Purpose

Amplitude envelopes of pulse waveforms. Syntax: waveform=pulse_shape(pulse_name,npoints)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- pulse_name -the name of the pulse (see function text)
- npoints -number of points in the pulse
- Output:
- waveform -normalised waveform as a vector

## Implementation structure

- Amplitude envelopes of pulse waveforms. Syntax:
- waveform=pulse_shape(pulse_name,npoints)
- pulse_name -the name of the pulse (see function text)
- npoints -number of points in the pulse
- Output:
- waveform -normalised waveform as a vector
- Check consistency
- Choose the shape
- Complain and bomb out
- Consistency enforcement
- Let me start with a parable. It concerns an Eastern European country
- whose parliament was considering a total smoking ban. In response, a

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `normpdf()`, `sinc()`, `ischar()`.
