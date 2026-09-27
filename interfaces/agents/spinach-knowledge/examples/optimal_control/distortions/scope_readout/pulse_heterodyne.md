# examples/optimal_control/distortions/scope_readout/pulse_heterodyne.m

- Signature: `pulse_heterodyne()`

## Purpose

Processes a recording of the proton optimal-control pulse, acquired on a 1 GHz oscilloscope using the probe's ¹³C coil, by heterodyning it to the rotating frame and comparing it with the ideal pulse.

## Method

Loads the scope data and ideal pulse, heterodynes the measured signal at 500.1029395 MHz, resamples its quadrature components by a factor of 1/20,000, and applies a 2.338 rad phase correction. The measured amplitude is scaled to the ideal amplitude; amplitude, in-phase, and out-of-phase components are then plotted for comparison.
