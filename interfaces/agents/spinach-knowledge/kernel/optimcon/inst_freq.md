# kernel/optimcon/inst_freq.m

- Signature: `freq=inst_freq(signal,dt,npoints,poly_order,amp_tol)`

## Purpose

Instantaneous frequency trajectory from a complex time-domain signal by regularised phase differentiation. Syntax: freq=inst_freq(signal,dt,npoints,poly_order,amp_tol)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- Instantaneous-frequency calculations differentiate pulse phase with respect to time. Numerically this is delicate because wrapped phases create 2π jumps, so implementations normally need unwrapping and careful finite-difference treatment to avoid spurious spikes.
- In pulse design, instantaneous frequency is the bridge between phase-modulated and frequency-modulated representations of a waveform; it is central for adiabatic sweeps, chirps, and hardware interpretability.

## Numerical / algorithmic content

- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.

## Parameters / inputs

- signal -complex row or column vector
- with time-domain signal
- dt -time step duration between
- signal points (seconds)
- npoints -odd number of signal points in the
- local least-squares window
- poly_order -local polynomial order used for
- phase differentiation
- amp_tol -fractional amplitude tolerance
- relative to the maximum signal
- amplitude; zero tolerance still
- masks zero-magnitude points where
- the phase is undefined

## Outputs

- freq -instantaneous frequency trajec-
- tory (Hz), same size as signal
- and on the same time grid
- Note: signal phase is unwrapped first and then differentiated
- by Savitzky-Golay local least-squares polynomial fits.
- This regularises numerical phase noise before the deri-
- vative is taken. Output points are set to NaN when any
- point in the local differentiation stencil is below the
- amplitude tolerance.

## Implementation structure

- Instantaneous frequency trajectory from a complex time-domain
- signal by regularised phase differentiation. Syntax:
- freq=inst_freq(signal,dt,npoints,poly_order,amp_tol)
- signal -complex row or column vector
- with time-domain signal
- dt -time step duration between
- signal points (seconds)
- npoints -odd number of signal points in the
- local least-squares window
- poly_order -local polynomial order used for
- phase differentiation
- amp_tol -fractional amplitude tolerance
