# experiments/esr_hyperfine/hyscore.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/esr_hyperfine/hyscore.m`
- Signature: `fid=hyscore(spin_system,parameters,H,R,K)`
- Total lines: 133

## Purpose

HYSCORE experiment, implemented as described in Szosenfogel and Goldfarb (http://dx.doi.org/10.1080/00268979809483260). Syntax: fid=hyscore(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Hyperfine ESR experiment implementations. These sequences probe coupled electron-nuclear dynamics through ENDOR or HYSCORE-type manipulations of coherence pathways.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters.nsteps number of points to be computed
- in each dimension
- parameters.sweep sweep width, Hz
- parameters.tau tau delay, seconds
- parameters.rho0 initial state
- parameters.coil detection state
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -two-dimensional free induction decay that Fourier
- transforms into a HYSCORE spectrum
- Note: the sequence uses ideal pulses, replace with shaped_pulse_af()
- to have soft pulses instead.

## Implementation structure

- HYSCORE experiment, implemented as described in Szosenfogel and
- Goldfarb (http://dx.doi.org/10.1080/00268979809483260). Syntax:
- fid=hyscore(spin_system,parameters,H,R,K)
- parameters.nsteps number of points to be computed
- in each dimension
- parameters.sweep sweep width, Hz
- parameters.tau tau delay, seconds
- parameters.rho0 initial state
- parameters.coil detection state
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `step()`, `evolution()`, `coherence()`, `ismatrix()`, `all()`, `isfield()`.
