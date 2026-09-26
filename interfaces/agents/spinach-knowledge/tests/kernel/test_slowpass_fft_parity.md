# tests/kernel/test_slowpass_fft_parity.m

- Signature: `result=test_slowpass_fft_parity()`

## Purpose

Tests slowpass amplitude normalisation against time-domain FFT. Syntax: result=test_slowpass_fft_parity()

## Physical / mathematical content

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Outputs

- result -regression test result with explanatory messages
- The test compares frequency-domain acquisition from slowpass() with the
- same damped one-spin signal acquired in the time domain and processed by
- Matlab's unnormalised FFT.

## Implementation structure

- Tests slowpass amplitude normalisation against time-domain FFT. Syntax:
- result=test_slowpass_fft_parity()
- result -regression test result with explanatory messages
- The test compares frequency-domain acquisition from slowpass() with the
- same damped one-spin signal acquired in the time domain and processed by
- Matlab's unnormalised FFT.
- Announce the test target
- State the slowpass normalisation target of the test
- Build a damped one-spin Liouville-space system
- Get production generators and states
- Acquire the signal in the time domain
- Use the exact FFT bins as the slowpass frequency grid
