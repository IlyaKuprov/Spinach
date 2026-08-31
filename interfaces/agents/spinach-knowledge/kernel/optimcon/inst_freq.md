# kernel/optimcon/inst_freq.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/inst_freq.m`
- Signature: `freq=inst_freq(signal,dt,npoints,poly_order,amp_tol)`
- Total lines: 126

## Purpose

Instantaneous frequency trajectory from a complex time-domain signal by regularised phase differentiation. Syntax: freq=inst_freq(signal,dt,npoints,poly_order,amp_tol)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- Instantaneous-frequency calculations differentiate pulse phase with respect to time. Numerically this is delicate because wrapped phases create 2π jumps, so implementations normally need unwrapping and careful finite-difference treatment to avoid spurious spikes.
- In pulse design, instantaneous frequency is the bridge between phase-modulated and frequency-modulated representations of a waveform; it is central for adiabatic sweeps, chirps, and hardware interpretability.

## Numerical / algorithmic content

- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 45-46: Check consistency; implemented by `grumble(signal,dt,npoints,poly_order,amp_tol)`.
- Lines 48-49: Input into column; implemented by `signal_col=signal(:)`.
- Lines 51-52: Unwrap the phase trajectory; implemented by `phase_col=unwrap(angle(signal_col))`.
- Lines 54-55: Differentiate the phase using local least-squares fits; implemented by `freq=sgolaydiff(phase_col,1,npoints,poly_order)/(2*pi*dt)`.
- Lines 57-58: Set the absolute amplitude threshold; implemented by `amp_limit=amp_tol*max(abs(signal_col))`.
- Lines 60-61: Preallocate the unreliable stencil mask; implemented by `weak_stencil=false(size(signal_col))`.
- Lines 63-64: Mark stencils that contain weak signal points; implemented by `win_half=(npoints-1)/2; nsamps=numel(signal_col)`.
- Lines 67-68: Select the differentiation stencil; implemented by `win_left=max(1,min(n-win_half,nsamps-npoints+1))`.
- Lines 71-72: Check whether any stencil point is below threshold; implemented by `weak_stencil(n)=any(abs(signal_col(win_idx))<=amp_limit)`.
- Lines 76-77: Remove unreliable instantaneous frequency values; implemented by `freq(weak_stencil)=NaN`.
- Lines 79-80: Shape back to input shape; implemented by `freq=reshape(freq,size(signal))`.

### Control flow inferred from the code

- Line 65: `for` loop over `n=1:nsamps`.

### Key state/data transformations

- Lines 49: computes `signal_col` using `signal_col=signal(:)`.
- Lines 52: computes `phase_col` using `phase_col=unwrap(angle(signal_col))`.
- Lines 55: computes `freq` using `freq=sgolaydiff(phase_col,1,npoints,poly_order)/(2*pi*dt)`.
- Lines 58: computes `amp_limit` using `amp_limit=amp_tol*max(abs(signal_col))`.
- Lines 61: computes `weak_stencil` using `weak_stencil=false(size(signal_col))`.
- Lines 64: computes `win_half` using `win_half=(npoints-1)/2; nsamps=numel(signal_col)`.
- Lines 68: computes `win_left` using `win_left=max(1,min(n-win_half,nsamps-npoints+1))`.
- Lines 69: computes `win_idx` using `win_idx=win_left:(win_left+npoints-1)`.
- Lines 77: computes `freq(weak_stencil)` using `freq(weak_stencil)=NaN`.

### Local helper functions

- Line 85: `grumble()` — `function grumble(signal,dt,npoints,poly_order,amp_tol)`.
  - Representative operation: `if (~isnumeric(signal))||(~isvector(signal))||isempty(signal)||isreal(signal)`.
  - Representative operation: `error('signal must be a non-empty complex numeric vector.')`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `signal()`, `unwrap()`, `sgolaydiff()`, `false()`, `weak_stencil()`, `any()`, `signal_col()`, `freq()`, `isvector()`, `isscalar()`.
