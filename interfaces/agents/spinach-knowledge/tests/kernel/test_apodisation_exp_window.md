# tests/kernel/test_apodisation_exp_window.m

- Signature: `result=test_apodisation_exp_window()`

## Purpose

Tests exponential FID apodisation. Syntax: result=test_apodisation_exp_window()

## Physical / mathematical content

- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test checks the explicit exponential window and the NMR Fourier
- convention that halves the first point of each active FID dimension.

## Implementation structure

- Tests exponential FID apodisation. Syntax:
- result=test_apodisation_exp_window()
- result -regression test result with explanatory messages
- The test checks the explicit exponential window and the NMR Fourier
- convention that halves the first point of each active FID dimension.
- Announce the test target
- State the processing target of the test
- Build a minimal reporting object and a constant FID
- Apply an exponential window
- Check the explicit window
