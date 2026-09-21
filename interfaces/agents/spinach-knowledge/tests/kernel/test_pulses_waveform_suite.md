# tests/kernel/test_pulses_waveform_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_pulses_waveform_suite.m`
- Signature: `result=test_pulses_waveform_suite()`
- Total lines: 126

## Purpose

Tests deterministic pulse waveform generators. Syntax: result=test_pulses_waveform_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test checks analytic waveform formulas, periodic phase tables,
- JCAMP pulse-file reading, chirp construction, and sech pulse coordinates.

## Implementation structure

- Tests deterministic pulse waveform generators. Syntax:
- result=test_pulses_waveform_suite()
- result -regression test result with explanatory messages
- The test checks analytic waveform formulas, periodic phase tables,
- JCAMP pulse-file reading, chirp construction, and sech pulse coordinates.
- Announce the test target
- State the waveform target of the test
- Check sawtooth and triangular wave formulae at simple fractions of the period
- Check Uhrig delay formula for three pulses over a unit interval
- Check periodic phase tables and wrap-around indexing
- Check simple analytic pulse envelopes
- Check reading of a distributed rectangular Bruker pulse file

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_close()`, `sawtooth()`, `triwave()`, `diff()`, `uhrig_times()`, `pmlg5()`, `spinal()`, `pulse_shape()`, `read_wave()`, `vg_pulse()`, `chirp_pulse()`, `polar2cartesian()`, `sech_pulse()`, `sech()`, `cosh()`.
