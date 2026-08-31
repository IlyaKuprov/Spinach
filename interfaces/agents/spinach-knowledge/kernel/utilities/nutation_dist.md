# kernel/utilities/nutation_dist.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/nutation_dist.m`
- Signature: `[freq,distr]=nutation_dist(curve,dt,lambda)`
- Total lines: 225

## Purpose

Nutation frequency distribution from a nutation curve measured with the same coil used for excitation and detection. Syntax: [freq,distr]=nutation_dist(curve,dt,lambda)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `fit_nutation()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 57-58: Check consistency; implemented by `grumble(curve,dt,lambda)`.
- Lines 60-61: Fold the input into a column; implemented by `signal=curve(:)`.
- Lines 65-66: Normalise the receiver gain away; implemented by `signal=signal/norm(signal,inf)`.
- Lines 68-69: Estimate the complex noise variance from the second difference; implemented by `diff_two=signal(3:end)-2*signal(2:end-1)+signal(1:end-2)`.
- Lines 75-76: Rotate the receiver phase line onto the real axis; implemented by `line_phase=exp(1i*angle(sum(signal.^2))/2)`.
- Lines 79-80: Apply a fade-out window for support estimation; implemented by `window=0.5+0.5*cos(pi*(0:(npts-1)).'/(npts-1))`.
- Lines 82-83: Sine-transform power spectrum on a zero-filled grid; implemented by `nfft=2^nextpow2(8*npts)`.
- Lines 90-91: Locate the noise-limited spectral support; implemented by `freq_step=2*pi/(nfft*dt)`.
- Lines 109-110: Double the support about its centre to expose the distribution margins; implemented by `freq_mid=(freq_lo+freq_hi)/2`.
- Lines 115-116: Choose a frequency grid at approximately twice the Fourier resolution; implemented by `resolution=ceil((freq_hi-freq_lo)*npts*dt/(2*pi))`.
- Lines 120-121: Second-difference regularisation matrix with clamped boundaries; implemented by `D=spdiags(ones(ngrid,1)*[1 -2 1],-1:1,ngrid,ngrid)`.
- Lines 123-124: Dimensionless reciprocity reception weight; implemented by `rec_weight=freq.'/freq_hi`.
- Lines 126-127: Silence the non-negative least squares solver; implemented by `options=optimset('Display','off')`.
- Lines 129-130: Search for a small error in the time origin; implemented by `shift_grid=linspace(-2*dt,2*dt,9)`.
- Lines 142-143: Refine the best time origin on the local interval; implemented by `shift_step=shift_grid(2)-shift_grid(1)`.
- Lines 154-155: Final fit at the corrected time origin; implemented by `kernel=sin((time+best_shift)*freq.').*rec_weight`.
- Lines 158-159: Convert fitted probability masses into a unit-integral density; implemented by `distr=max(best_weights,0)`.

### Control flow inferred from the code

- Line 95: conditional branch on `isempty(active)`.
- Line 101: conditional branch on `freq_hi<=freq_lo`.
- Line 105: conditional branch on `freq_hi<=freq_lo`.
- Line 133: `for` loop over `n=1:numel(shift_grid)`.
- Line 136: conditional branch on `fit_err<best_err`.
- Line 145: `for` loop over `n=1:numel(shift_grid)`.
- Line 148: conditional branch on `fit_err<best_err`.

### Key state/data transformations

- Lines 61: computes `signal` using `signal=curve(:)`.
- Lines 62: computes `npts` using `npts=numel(signal)`.
- Lines 63: computes `time` using `time=(0:(npts-1)).'*dt`.
- Lines 69: computes `diff_two` using `diff_two=signal(3:end)-2*signal(2:end-1)+signal(1:end-2)`.
- Lines 70: computes `noise_var` using `noise_var=median(abs(diff_two).^2)/(6*log(2))`.
- Lines 71: computes `signal_peak` using `signal_peak=max(abs(signal).^2)`.
- Lines 76: computes `line_phase` using `line_phase=exp(1i*angle(sum(signal.^2))/2)`.
- Lines 77: computes `rotated` using `rotated=real(conj(line_phase)*signal)`.
- Lines 80: computes `window` using `window=0.5+0.5*cos(pi*(0:(npts-1)).'/(npts-1))`.
- Lines 83: computes `nfft` using `nfft=2^nextpow2(8*npts)`.
- Lines 84: computes `spectrum` using `spectrum=fftshift(fft(window.*rotated,nfft))`.
- Lines 85: computes `freq_axis` using `freq_axis=(-floor(nfft/2):ceil(nfft/2)-1).'*2*pi/(nfft*dt)`.
- Lines 87: computes `selected_freq` using `selected_freq=freq_axis(selected)`.
- Lines 88: computes `selected_power` using `selected_power=imag(spectrum(selected)).^2`.
- Lines 91: computes `freq_step` using `freq_step=2*pi/(nfft*dt)`.
- Lines 92: computes `noise_bin` using `noise_bin=noise_var*sum(window.^2)/4`.
- Lines 93: computes `power_cut` using `power_cut=max(10*noise_bin,1e-5*max(selected_power))`.
- Lines 97: computes `active` using `active=max(1,peak_idx-2):min(numel(selected_freq),peak_idx+2)`.

### Local helper functions

- Line 165: `fit_nutation()` — `function [weights,fit_err]=fit_nutation(signal,kernel,D,lambda,options)`. Stack the data fit and the regularisation penalty
  - Representative operation: `lsq_mat=[kernel;sqrt(lambda)*full(D)]`.
  - Representative operation: `zero_pen=zeros(size(D,1),1)`.
- Line 203: `grumble()` — `function grumble(curve,dt,lambda)`.
  - Representative operation: `if (~isnumeric(curve))||(~isvector(curve))||isempty(curve)`.
  - Representative operation: `error('curve must be a non-empty numeric vector.')`.

## Parameters / inputs

- curve -nutation curve, a row or column vector; either
- the complex X+iY output of a quadrature receiver,
- or a real phase-corrected trace
- dt -sampling interval in seconds
- lambda -second-derivative Tikhonov regularisation para-
- meter, a non-negative real scalar; the curve is
- normalised to unit maximum modulus and the fit-
- ting kernel is dimensionless, so lambda is a di-
- mensionless number of the order of the ratio of
- the squared Frobenius norms of the kernel and of
- the second difference matrix; zero switches the
- regularisation off

## Outputs

- freq -nutation frequency grid in rad/s, a column vector
- distr -non-negative nutation frequency density in inverse
- rad/s, normalised to unit integral over freq
- Notes: when one coil both excites and detects, the reciprocity
- principle makes the detected amplitude of every isochromat
- proportional to its own nutation frequency, and only the
- sine component of the nutation is observable. The curve is
- therefore modelled as an unknown complex receiver scale
- times
- s(t)=integral(distr(freq)*freq*sin(freq*(t+t0)),d freq)
- and the reception weight is divided out, so that the
- returned density is the true nutation frequency distribu-
- tion. A real receiver scale is a valid special case, and
- a phase-corrected real trace is therefore accepted. The
- frequency support, receiver phase, and sub-sample time
- shift are selected from the supplied trace; the returned
- density is therefore a stable estimate rather than a raw
- finite-time transform. The reconstruction runs on twice
- the noise-limited support width, centred on the same band
- and clipped to the Nyquist interval, so that the margins
- of the distribution are visible rather than truncated.

## Implementation structure

- Nutation frequency distribution from a nutation curve measured with
- the same coil used for excitation and detection. Syntax:
- [freq,distr]=nutation_dist(curve,dt,lambda)
- curve -nutation curve, a row or column vector; either
- the complex X+iY output of a quadrature receiver,
- or a real phase-corrected trace
- dt -sampling interval in seconds
- lambda -second-derivative Tikhonov regularisation para-
- meter, a non-negative real scalar; the curve is
- normalised to unit maximum modulus and the fit-
- ting kernel is dimensionless, so lambda is a di-
- mensionless number of the order of the ratio of

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `curve()`, `signal()`, `median()`, `conj()`, `nextpow2()`, `fftshift()`, `freq_axis()`, `spectrum()`, `selected_freq()`, `active()`, `spdiags()`, `optimset()`, `shift_grid()`, `fit_nutation()`, `trapz()`.
