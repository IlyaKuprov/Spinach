# kernel/derivatives/sgolaydiff.m

- Signature: `dy=sgolaydiff(y,der_order,npoints,poly_order)`

## Purpose

Savitzky-Golay differentiation of noisy sampled signals by local least-squares polynomial fitting.

## Physical / mathematical content

- Fits a local polynomial around each sample and evaluates the requested derivative at that sample. Derivative order 0 returns the smoothed signal.
- Derivatives are reported for a uniform grid with unit sample spacing.

## Numerical / algorithmic content

- Uses an odd-length window, with sided windows near the signal boundaries and centred windows elsewhere.
- Centres and scales the sample offsets, constructs a local Vandermonde matrix, and fits the polynomial by least squares using MATLAB's backslash operator.
- Multiplies the fitted coefficient by the derivative-order factorial and divides by the offset scale raised to that order.

## Parameters / inputs

- `y` — N-by-M signal matrix; rows are samples and columns are independent signals. It must be a non-empty, dense floating-point matrix with finite values and at least three rows.
- `der_order` — non-negative integer derivative order; must not exceed `poly_order`. Order 0 returns the smoothed signal.
- `npoints` — odd number of points in the local least-squares window; at least 3 and no greater than the number of sample rows.
- `poly_order` — non-negative integer order of the local polynomial; must be smaller than `npoints`.

## Outputs

- `dy` — N-by-M derivative matrix on a unit-step uniform grid.

## Note

- `sgolaydiff(s,1,7,3)` is recommended for differentiating EPR spectra; use a tight integration tolerance and increase the number of field/frequency axis points.

## Reference

- https://spindynamics.org/wiki/index.php?title=sgolaydiff.m