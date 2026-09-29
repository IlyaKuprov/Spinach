# examples/optimal_control/pulse_analysis.m

- Signature: `pulse_analysis()`
- Source: [`examples/optimal_control/pulse_analysis.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/pulse_analysis.m)

## Purpose

A short signal-processing illustration of a quadratic-chirp superposition and its spectrogram; the source says it was adapted from the MATLAB example set and gives a calculation time of seconds. It does not define a spin system or perform a GRAPE pulse optimisation.

## Signal and analysis settings

At a sampling frequency of 1000 Hz, the script samples two seconds from 0 through 1.999 s. It adds two quadratic chirps: one starts at 100 Hz and reaches 200 Hz at 1 s, while the other starts at 200 Hz and reaches 100 Hz at 1 s. The first subplot shows signal amplitude in arbitrary units versus time. The second displays a spectrogram using a 100-sample window, 80-sample overlap, 100 frequency points, and a minimum threshold of −50 dB; its axes are time in seconds and frequency in Hz.
