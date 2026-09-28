# kernel/utilities/nutation_dist.m

- Signature: `[freq,distr]=nutation_dist(curve,dt,lambda)`

## Purpose

Estimate a nutation frequency distribution from a nutation curve measured with the same coil for excitation and detection. The result is a non-negative, unit-integral density obtained by fitting the curve rather than taking a raw finite-time transform.

## Inputs

- `curve` — nutation curve as a row or column vector: either the complex `X+iY` output of a quadrature receiver or a real phase-corrected trace. It must be numeric, finite, nonzero, and contain at least eight points.
- `dt` — sampling interval in seconds; a positive real scalar.
- `lambda` — second-derivative Tikhonov regularisation parameter; a non-negative real scalar. Zero disables regularisation. The curve is normalised to unit maximum modulus and the fitting kernel is dimensionless, so `lambda` is dimensionless and is of the order of the ratio of the squared Frobenius norms of the kernel and the second-difference matrix.

## Outputs

- `freq` — nutation frequency grid in rad/s, returned as a column vector.
- `distr` — non-negative nutation frequency density in inverse rad/s, normalised to unit integral over `freq`.

## Model and method

When one coil both excites and detects, reciprocity makes the detected amplitude of each isochromat proportional to its nutation frequency; only the sine component is observable. The measured curve is modelled as an unknown complex receiver scale times `s(t)=integral(distr(freq)*freq*sin(freq*(t+t0)),d freq)`. The fit divides out this reception weight to recover the nutation frequency distribution. A real receiver scale is a valid special case, allowing a phase-corrected real trace.

The routine normalises the signal, estimates noise from its second difference, and uses a faded, zero-filled Fourier spectrum to identify noise-limited frequency support. It doubles that support width about the same centre, subject to the `0` to `pi/dt` Nyquist interval, so the distribution margins remain visible. A frequency grid of 160 to 400 points is then fitted with non-negative least squares and a second-difference Tikhonov penalty. The fit selects receiver phase and a sub-sample time-origin shift from the supplied trace before converting fitted weights to a unit-integral density. If the curve contains no resolvable frequency range, the routine reports an error.

Source: <https://spindynamics.org/wiki/index.php?title=nutation_dist.m>