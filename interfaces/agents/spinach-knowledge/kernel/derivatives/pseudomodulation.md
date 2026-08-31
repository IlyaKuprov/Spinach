# kernel/derivatives/pseudomodulation.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/derivatives/pseudomodulation.m`
- Signature: `output=pseudomodulation(field,spectrum,mod_amp,mod_order)`
- Total lines: 120

## Purpose

Pseudomodulation of uniformly sampled spectra using the Hyde et al. Fourier-domain algorithm. Syntax: output=pseudomodulation(field,spectrum,mod_amp,mod_order)

## Physical / mathematical content

- Derivative utilities. These routines compute finite-difference, analytical, or optimisation-oriented derivatives needed for sensitivity analysis, fitting, and optimal control.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 34-35: Check consistency; implemented by `grumble(field,spectrum,mod_amp,mod_order)`.
- Lines 37-38: Get field grid parameters; implemented by `npts=numel(field)`.
- Lines 41-42: Build the angular frequency axis in Matlab FFT ordering; implemented by `ang_freq=sign(field_step)*real(fftdiff(1,npts,abs(field_step))/1i)`.
- Lines 45-46: Convert the modulation amplitude into Bessel arguments; implemented by `bessel_arg=mod_amp*ang_freq/2`.
- Lines 48-49: Fourier transform along the field axis; implemented by `spec_ft=fft(spectrum,npts,1)`.
- Lines 51-52: Apply the requested pseudomodulation harmonic; implemented by `switch mod_order`.
- Lines 56-57: Zeroth harmonic; implemented by `output=ifft(spec_ft.*besselj(0,bessel_arg),npts,1)`.
- Lines 61-62: First harmonic; implemented by `output=2i*ifft(spec_ft.*besselj(1,bessel_arg),npts,1)`.
- Lines 66-67: Second harmonic; implemented by `output=2*ifft(spec_ft.*besselj(2,bessel_arg),npts,1)`.
- Lines 71-72: Remove round-off imaginary parts from real input spectra; implemented by `if isreal(spectrum)`.

### Control flow inferred from the code

- Line 52: dispatches on `mod_order`; cases `0`, `1`, `2`.
- Line 72: conditional branch on `isreal(spectrum)`.

### Key state/data transformations

- Lines 38: computes `npts` using `npts=numel(field)`.
- Lines 39: computes `field_step` using `field_step=field(2)-field(1)`.
- Lines 42: computes `ang_freq` using `ang_freq=sign(field_step)*real(fftdiff(1,npts,abs(field_step))/1i)`.
- Lines 46: computes `bessel_arg` using `bessel_arg=mod_amp*ang_freq/2`.
- Lines 49: computes `spec_ft` using `spec_ft=fft(spectrum,npts,1)`.
- Lines 57: computes `output` using `output=ifft(spec_ft.*besselj(0,bessel_arg),npts,1)`.

### Local helper functions

- Line 79: `grumble()` — `function grumble(field,spectrum,mod_amp,mod_order)`.
  - Representative operation: `if (~isfloat(field))||isempty(field)||issparse(field)||(~iscolumn(field))`.
  - Representative operation: `error('field must be a non-empty dense floating-point column vector.')`.

## Parameters / inputs

- field -N-by-1 real, ordered, uniformly spaced field
- axis
- spectrum -N-by-M spectrum matrix; rows are field samples,
- and columns are independent spectra
- mod_amp -non-negative modulation amplitude in field units
- mod_order -modulation harmonic order: 0, 1, or 2

## Outputs

- output -N-by-M pseudomodulated spectrum matrix
- The implementation follows Eqs. 5-7 of Hyde et al., J. Magn.
- Reson. 96, 1-13 (1992). After phase-sensitive detection, the
- time-dependent prefactors are set to unity, leaving amplitude
- factors 2i for the first harmonic, and 2 for the second harmonic.

## Implementation structure

- Pseudomodulation of uniformly sampled spectra using the Hyde
- et al. Fourier-domain algorithm. Syntax:
- output=pseudomodulation(field,spectrum,mod_amp,mod_order)
- field -N-by-1 real, ordered, uniformly spaced field
- axis
- spectrum -N-by-M spectrum matrix; rows are field samples,
- and columns are independent spectra
- mod_amp -non-negative modulation amplitude in field units
- mod_order -modulation harmonic order: 0, 1, or 2
- output -N-by-M pseudomodulated spectrum matrix
- The implementation follows Eqs. 5-7 of Hyde et al., J. Magn.
- Reson. 96, 1-13 (1992). After phase-sensitive detection, the

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `field()`, `sign()`, `fftdiff()`, `ang_freq()`, `besselj()`, `isfloat()`, `issparse()`, `iscolumn()`, `any()`, `eps()`, `ismatrix()`, `spectrum()`, `ismember()`.
