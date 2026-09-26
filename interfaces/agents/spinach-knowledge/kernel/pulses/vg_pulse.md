# kernel/pulses/vg_pulse.m

- Signature: `waveform=vg_pulse(pulse_name,npoints,duration)`

## Purpose

Veshtort-Griffin shaped pulses, generated from tables given in There are good reasons to believe (see Section 2.2 of the paper) that these are the best possible pulses within their design specifications and basis sets. Syntax: waveform=vg_pulse(pulse_name,npoints,duration)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

## Parameters / inputs

- pulse_name -a character string, one of the following: E0A,
- E0B, E100A, E100B, E200A, E200D, E200F, E300C,
- E300F, E400B, E300A, E500A, E500B, E500C, E600A,
- E600C, E600F, E800A, E800B, E1000B
- npoints -number of discrete time intervals in the pulse
- duration -duration of the pulse, seconds

## Outputs

- waveform -amplitude of the pulse at each interval (there
- is no phase modulation), normalised to produce
- a 90-degree pulse, rad/s

## Implementation structure

- Veshtort-Griffin shaped pulses, generated from tables given in
- There are good reasons to believe (see Section 2.2 of the paper) that
- these are the best possible pulses within their design specifications
- and basis sets. Syntax:
- waveform=vg_pulse(pulse_name,npoints,duration)
- pulse_name -a character string, one of the following: E0A,
- E0B, E100A, E100B, E200A, E200D, E200F, E300C,
- E300F, E400B, E300A, E500A, E500B, E500C, E600A,
- E600C, E600F, E800A, E800B, E1000B
- npoints -number of discrete time intervals in the pulse
- duration -duration of the pulse, seconds
- waveform -amplitude of the pulse at each interval (there
