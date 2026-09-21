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
