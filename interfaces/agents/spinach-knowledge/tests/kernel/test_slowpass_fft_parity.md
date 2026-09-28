# tests/kernel/test_slowpass_fft_parity.m

- Signature: `result=test_slowpass_fft_parity()`

## Purpose

Checks slowpass amplitude normalisation against a time-domain signal processed with MATLAB's unnormalised FFT.

## Physical / mathematical content

For a damped one-spin signal, frequency-domain acquisition by `slowpass()` should match the unnormalised FFT amplitude at zero frequency.

## Numerical / algorithmic content

The test builds a damped one-proton Liouville-space system at `14.1 T`, acquires a `4096`-point FID with a `4096 Hz` sweep, and compares the zero-frequency bin (index `2049`) of `slowpass()` with `fftshift(fft(fid))`. The comparison uses absolute and relative tolerances of `1e-6` and `2e-3`. It also verifies that a one-point spectrum is rejected when the error message contains `integer greater than one`.

## Outputs

`result` is the regression-test record with the comparison and rejection checks.

## Implementation structure

Acquires the same damped signal in time and frequency domains, checks the zero-frequency amplitude, and tests the single-point error path.

## Header notes

The source header credits Ilya Kuprov (`ilya.kuprov@weizmann.ac.il`).
