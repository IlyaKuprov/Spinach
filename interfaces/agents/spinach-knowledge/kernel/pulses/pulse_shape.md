# kernel/pulses/pulse_shape.m

- Signature: `waveform=pulse_shape(pulse_name,npoints)`

## Purpose

Amplitude envelopes of pulse waveforms.


## Parameters / inputs

- pulse_name -the name of the pulse (see function text)
- npoints -number of points in the pulse
- Output:
- waveform -normalised waveform as a vector

## Implementation structure

- Check consistency.
- Choose the pulse shape.
- Reject invalid pulse shapes.
- Enforce consistency of the output.
