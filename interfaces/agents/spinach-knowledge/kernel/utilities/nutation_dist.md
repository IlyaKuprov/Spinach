# kernel/utilities/nutation_dist.m

## Purpose

Estimates the nutation frequency distribution from a nutation curve acquired with the same radiofrequency coil used for both excitation and detection, as documented in the file header. The function returns a non-negative frequency density rather than a raw finite-time transform.

## Behaviour

- Syntax: `[freq,distr]=nutation_dist(curve,dt,lambda)`.
- Input validation (`grumble`) requires `curve` to be a non-empty numeric vector with at least eight points, no `Inf` or `NaN` entries, and not identically zero; `dt` must be a positive real scalar; `lambda` must be a non-negative real scalar.
- The curve is folded into a column, a time axis is built as `(0:(npts-1)).'*dt`, and the signal is normalised by its infinity norm.
- Noise variance is estimated from the second difference as `median(abs(diff_two).^2)/(6*log(2))`, clipped to at most `0.05` of the signal peak power and at least `eps*norm(signal,2)^2/npts`.
- The receiver phase line is rotated onto the real axis using `exp(1i*angle(sum(signal.^2))/2)`.
- A cosine fade-out window `0.5+0.5*cos(pi*(0:(npts-1)).'/(npts-1))` is applied before a zero-filled FFT on a grid of size `2^nextpow2(8*npts)`; the non-negative-frequency sine-transform power `imag(spectrum).^2` is used to locate the noise-limited support, with a power cut of `max(10*noise_bin,1e-5*max(selected_power))` where `noise_bin=noise_var*sum(window.^2)/4`.
- The support is extended by two frequency steps (`2*pi/(nfft*dt)`) on each side, then doubled about its centre and clipped to `[0, pi/dt]` so the distribution margins remain visible; if no resolvable range exists the function errors with `'the curve does not contain a resolvable frequency range.'`.
- The frequency grid uses `ngrid=min(400,max(160,2*resolution))` points, where `resolution=ceil((freq_hi-freq_lo)*npts*dt/(2*pi))`.
- A second-difference regularisation matrix `D=spdiags(ones(ngrid,1)*[1 -2 1],-1:1,ngrid,ngrid)` is used with the Tikhonov parameter `lambda`.
- The reciprocity reception weight `freq.'/freq_hi` multiplies the sine kernel `sin((time+shift)*freq.')`, and the reception weight is divided out of the returned density.
- A sub-sample time-origin search runs over `linspace(-2*dt,2*dt,9)` and then refines over a five-point local interval of width one coarse step; the shift with the lowest fit error is used for the final fit.
- `fit_nutation` stacks `[kernel;sqrt(lambda)*full(D)]` with a zero penalty vector, tries both signs of the receiver phase estimate, alternates phase estimation with `lsqnonneg` for four iterations per branch, and keeps the branch with the smaller squared residual `norm(phase*predicted-signal,2)^2`.
- The final weights are clamped to non-negative values and normalised with `trapz` so that `distr` integrates to one over `freq`.

## Inputs and outputs

**Inputs**

- `curve` — nutation curve, a row or column vector; either the complex `X+iY` output of a quadrature receiver or a real phase-corrected trace.
- `dt` — sampling interval in seconds.
- `lambda` — second-derivative Tikhonov regularisation parameter, a non-negative real scalar. Because the curve is normalised to unit maximum modulus and the fitting kernel is dimensionless, `lambda` is a dimensionless number of order the ratio of the squared Frobenius norms of the kernel and the second-difference matrix; zero switches regularisation off.

**Outputs**

- `freq` — nutation frequency grid in rad/s, a column vector.
- `distr` — non-negative nutation frequency density in inverse rad/s, normalised to unit integral over `freq`.

## References

- Source: [kernel/utilities/nutation_dist.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/nutation_dist.m)
- Spinach Wiki: [nutation_dist.m](https://spindynamics.org/wiki/index.php?title=nutation_dist.m)
