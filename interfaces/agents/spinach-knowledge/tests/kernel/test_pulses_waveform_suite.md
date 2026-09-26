# tests/kernel/test_pulses_waveform_suite.m

- Signature: `result=test_pulses_waveform_suite()`

## Purpose

Tests deterministic pulse waveform generators. Syntax: result=test_pulses_waveform_suite()

## Physical / mathematical content

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
