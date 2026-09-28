# kernel/optimcon/distortions/firf.m

- Signature: `[w,J]=firf(w,ker)`

## Purpose

Applies an FIR convolution filter to a Spinach optimal-control waveform. Each pair of rows represents the in-phase and quadrature components of one complex control channel. The distal end of the convolution is truncated to retain the input number of time samples.

## Parameters / inputs

- `w`: Real numeric waveform with one time slice per column and rows arranged `XYXY...` across control channels. The number of rows must be even.
- `ker`: Nonempty numeric vector of FIR filter coefficients. Coefficients may be complex.

## Outputs

- `w`: Filtered waveform with the same dimensions as the input. Leaving sufficient ring-down margin is the user's responsibility.
- `J`: When requested, the sparse Jacobian with respect to vectorisations of the output and input arrays.

## Numerical / algorithmic content

The function constructs a sparse Toeplitz convolution matrix, using the available coefficients up to the waveform length. For each channel, it combines the paired rows into a complex signal, applies the filter, and writes the real and imaginary parts back to the corresponding rows. When requested, it assembles the Jacobian from the real and imaginary parts of the convolution matrix.

## Reference

- <https://spindynamics.org/wiki/index.php?title=firf.m>