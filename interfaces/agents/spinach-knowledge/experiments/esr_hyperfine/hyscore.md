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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 40-41: Check consistency; implemented by `grumble(parameters,H,R,K)`.
- Lines 43-44: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 46-47: Get the pulse operators; implemented by `Lp=operator(spin_system,'L+','E')`.
- Lines 51-52: Apply the first pulse; implemented by `rho=step(spin_system,Lx,parameters.rho0,pi/2)`.
- Lines 54-55: Run the tau evolution; implemented by `rho=evolution(spin_system,L,[],rho,parameters.tau,1,'final')`.
- Lines 57-58: Apply the second pulse; implemented by `rho=step(spin_system,Lx,rho,pi/2)`.
- Lines 60-61: Apply coherence filter; implemented by `rho=coherence(spin_system,rho,{{'E',0}})`.
- Lines 63-65: Run the indirect dimension evolution; implemented by `rho_stack=evolution(spin_system,L,[],rho,1/parameters.sweep, parameters.nsteps(1)-1,'trajectory')`.
- Lines 67-68: Apply the third pulse; implemented by `rho_stack=step(spin_system,Lx,rho_stack,pi)`.
- Lines 70-71: Propagate coil state backwards in time; implemented by `coil=evolution(spin_system,L,[],parameters.coil,-parameters.tau,1,'final')`.
- Lines 73-74: Apply a backwards pulse on the coil; implemented by `coil=step(spin_system,-Lx,coil,pi/2)`.
- Lines 76-78: Detect on new coil state in the direct dimension; implemented by `fid=evolution(spin_system,L,coil,rho_stack,1/parameters.sweep, parameters.nsteps(2)-1,'observable')`.

### Key state/data transformations

- Lines 44: computes `L` using `L=H+1i*R+1i*K`.
- Lines 47: computes `Lp` using `Lp=operator(spin_system,'L+','E')`.
- Lines 48: computes `Lm` using `Lm=operator(spin_system,'L-','E')`.
- Lines 49: computes `Lx` using `Lx=(Lp+Lm)/2`.
- Lines 52: computes `rho` using `rho=step(spin_system,Lx,parameters.rho0,pi/2)`.
- Lines 64-65: computes `rho_stack` using `rho_stack=evolution(spin_system,L,[],rho,1/parameters.sweep, parameters.nsteps(1)-1,'trajectory')`.
- Lines 71: computes `coil` using `coil=evolution(spin_system,L,[],parameters.coil,-parameters.tau,1,'final')`.
- Lines 77-78: computes `fid` using `fid=evolution(spin_system,L,coil,rho_stack,1/parameters.sweep, parameters.nsteps(2)-1,'observable')`.

### Local helper functions

- Line 83: `grumble()` — `function grumble(parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))|| (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.
  - Representative operation: `(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.

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
