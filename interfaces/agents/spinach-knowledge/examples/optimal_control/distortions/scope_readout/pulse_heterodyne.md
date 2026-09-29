# examples/optimal_control/distortions/scope_readout/pulse_heterodyne.m

Source: [examples/optimal_control/distortions/scope_readout/pulse_heterodyne.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/distortions/scope_readout/pulse_heterodyne.m)

- Signature: `pulse_heterodyne()`

## Purpose and inputs

This script processes an imported oscilloscope trace of a proton optimal-control pulse and compares it visually with the ideal pulse waveform. The source comment describes the recording as made with a 1-GHz oscilloscope through the same probe's 13C coil. It loads `time_grid` and `expt_data` from `scope_readout.mat`, and loads `pulse` from `ideal_pulse.mat`; it does not simulate an oscilloscope or generate the measured input. The imported trace is therefore distinct from the ideal control waveform and from the Spinach simulations in the RLC examples.

## Processing and plots

The script plots the raw trace against its wall-clock time (seconds), then uses the first time-grid interval as the sample spacing and heterodynes the trace at 500.1029395 MHz. It separately resamples the returned real and imaginary components by 1/20,000, rebuilds a time grid over the original endpoints, and forms the complex signal. An empirical phase correction of 2.338 radians is applied. The measured complex trace is scaled so its peak amplitude matches that of the ideal waveform; the ideal pulse's two rows are interpreted as the real and imaginary controls, and a zero column is appended before plotting on a 0-to-0.1-second theoretical grid.

The figure compares the envelope, in-phase component, and out-of-phase component, with the component plots limited to 0–0.05 seconds. These are visual comparisons after the stated phase and peak-amplitude adjustments; the script does not calculate a fit score, report calibrated field values, or itself establish hardware validation.
