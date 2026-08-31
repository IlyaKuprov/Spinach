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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 21-22: Check consistency; implemented by `grumble(pulse_name,npoints)`.
- Lines 24-25: Choose the shape; implemented by `switch pulse_name`.
- Lines 48-49: Complain and bomb out; implemented by `error('unknown pulse name.')`.

### Control flow inferred from the code

- Line 25: dispatches on `pulse_name`; cases `'gaussian'`, `'sinc5'`, `'sinc3'`, `'rectangular'`.

### Key state/data transformations

- Lines 29: computes `time_grid` using `time_grid=linspace(-2,2,npoints)`.
- Lines 30: computes `waveform` using `waveform=normpdf(time_grid)/sqrt(2)`.

### Local helper functions

- Line 56: `grumble()` — `function grumble(pulse_name,npoints)`. Let me start with a parable. It concerns an Eastern European country whose parliament was considering a total smoking ban. In response, a
  - Representative operation: `if ~ischar(pulse_name)`.
  - Representative operation: `error('pulse_name parameter must be a character string.')`.

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
