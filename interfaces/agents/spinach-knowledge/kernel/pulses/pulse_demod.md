# kernel/pulses/pulse_demod.m

- Signature: `demod_pulse=pulse_demod(time_grid,in_phase,out_phase)`

## Purpose

Interactive demodulation of a complex pulse waveform by a user- specified frequency. Syntax: demod_pulse=pulse_demod(time_grid,in_phase,out_phase)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

## Parameters / inputs

- time_grid -strictly increasing time grid, seconds
- in_phase -in-phase pulse waveform component
- out_phase -out-of-phase pulse waveform component

## Outputs

- demod_pulse -demodulated complex pulse waveform

## Header notes

- The frequency entry field uses Hz. The GHz, MHz, kHz, and Hz buttons
- set the slider step size without changing the demodulation frequency.
- The phase and frequency buttons switch the plot between unwrapped
- phase in radians and instantaneous frequency in Hz. The sticky wrap
- button switches the phase plot into the [0,2*pi] interval.
- The slider range is a moving 100-step window in the selected units.
- The save button returns the current demodulated waveform and exits.
- The complex waveform in_phase+1i*out_phase is multiplied by
- exp(2*pi*1i*freq*time_grid).

## Implementation structure

- Interactive demodulation of a complex pulse waveform by a user-
- specified frequency. Syntax:
- demod_pulse=pulse_demod(time_grid,in_phase,out_phase)
- time_grid -strictly increasing time grid, seconds
- in_phase -in-phase pulse waveform component
- out_phase -out-of-phase pulse waveform component
- demod_pulse -demodulated complex pulse waveform
- The frequency entry field uses Hz. The GHz, MHz, kHz, and Hz buttons
- set the slider step size without changing the demodulation frequency.
- The phase and frequency buttons switch the plot between unwrapped
- phase in radians and instantaneous frequency in Hz. The sticky wrap
- button switches the phase plot into the [0,2*pi] interval.
