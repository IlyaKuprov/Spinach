# examples/optimal_control/pulse_analysis.m

- Signature: `pulse_analysis()`

## Purpose

An example of spectrogram analysis for a quadratic chirp pulse; adapted from Matlab example set. Calculation time: seconds.

## Physical / mathematical content

- The signal is a superposition of two quadratic chirps: one sweeps from 100 to 200 Hz and the other from 200 to 100 Hz over one second.

## Numerical / algorithmic content

- Samples the signal at 1000 Hz over a two-second time vector.
- Computes a spectrogram with a 100-sample window, 80-sample overlap, 100 frequency points, and a −50 dB minimum threshold.

## Implementation structure

- Plots the signal amplitude versus time and its spectrogram versus time and frequency in two subplots.
