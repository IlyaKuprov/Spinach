# kernel/optimcon/distortions/firf.m

- Signature: [w,J]=firf(w,ker)
- MATLAB source: [firf.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/distortions/firf.m)

## Purpose

Applies a finite-impulse-response convolution filter to each complex control channel in a Spinach optimal-control waveform. The output keeps the input sample count by truncating the distal end of the convolution; leaving sufficient ring-down margin is the user's responsibility.

## Inputs and units

- w is a real numeric waveform. Columns are time samples; rows are arranged XYXY..., pairing each channel's in-phase X and quadrature Y components. The number of rows must be even.
- ker is a nonempty numeric vector of FIR coefficients. The implementation accepts complex coefficients as well as real ones. Coefficients multiply waveform samples; no separate coefficient units are specified by the function.

## Output and derivative

For each channel the paired rows are combined as z=X+iY, filtered by a Toeplitz convolution matrix, then split back into real and imaginary rows. The matrix is lower-triangular in its causal sample ordering. If ker has fewer entries than the waveform has time samples, it is zero-padded to that length; if it has more, only the first number-of-samples entries are used. The returned w has the same dimensions as the input.

The optional J is a sparse Jacobian of the vectorised output with respect to the vectorised input. It represents the real-coordinate derivative of the complex linear filter, including the X/Y cross terms when the filter coefficients are complex. No adjoint is returned.

## References

- Spinach documentation: https://spindynamics.org/wiki/index.php?title=firf.m
