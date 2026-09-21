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
