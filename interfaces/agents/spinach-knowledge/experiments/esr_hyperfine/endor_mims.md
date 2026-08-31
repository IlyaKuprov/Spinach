# experiments/esr_hyperfine/endor_mims.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/esr_hyperfine/endor_mims.m`
- Signature: `fid=endor_mims(spin_system,parameters,H,R,K)`
- Total lines: 114

## Purpose

Mims ENDOR pulse sequence with ideal hard pulses. Syntax: fid=endor_mims(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Hyperfine ESR experiment implementations. These sequences probe coupled electron-nuclear dynamics through ENDOR or HYSCORE-type manipulations of coherence pathways.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 30-31: Move into adjoint representation if needed; implemented by `[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)`.
- Lines 33-34: Consistency check; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 36-37: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 39-40: Initial state; implemented by `rho=state(spin_system,'Lz','electrons')`.
- Lines 42-43: Detection state; implemented by `coil=state(spin_system,'L+','electrons')`.
- Lines 45-46: Pulse operators; implemented by `Ep=operator(spin_system,'L+','electrons'); Ey=(Ep-Ep')/2i`.
- Lines 49-50: Apply the initial pulses; implemented by `rho=step(spin_system,Ey,rho,pi/2)`.
- Lines 54-55: Analytical coherence selection; implemented by `rho=coherence(spin_system,rho,{{'electrons',0}})`.
- Lines 57-58: Apply pulses on nuclear spins; implemented by `rho=step(spin_system,Ny,rho,pi/2)`.
- Lines 60-61: Analytical coherence selection; implemented by `rho=coherence(spin_system,rho,{{'nuclei',[-1 1]}})`.
- Lines 63-65: Indirect dimension evolution; implemented by `rho_stack=evolution(spin_system,L,[],rho,1/parameters.sweep, parameters.npoints-1,'trajectory')`.
- Lines 69-70: Apply pulse on electron spin to refocus; implemented by `rho_stack=step(spin_system,Ey,rho_stack,pi/2)`.
- Lines 74-75: Detect; implemented by `fid=full(transpose(coil'*rho_stack))`.

### Key state/data transformations

- Lines 31: computes `[spin_system,parameters,H,R,K]` using `[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)`.
- Lines 37: computes `L` using `L=H+1i*R+1i*K`.
- Lines 40: computes `rho` using `rho=state(spin_system,'Lz','electrons')`.
- Lines 43: computes `coil` using `coil=state(spin_system,'L+','electrons')`.
- Lines 46: computes `Ep` using `Ep=operator(spin_system,'L+','electrons'); Ey=(Ep-Ep')/2i`.
- Lines 47: computes `Np` using `Np=operator(spin_system,'L+','nuclei'); Ny=(Np-Np')/2i`.
- Lines 64-65: computes `rho_stack` using `rho_stack=evolution(spin_system,L,[],rho,1/parameters.sweep, parameters.npoints-1,'trajectory')`.
- Lines 75: computes `fid` using `fid=full(transpose(coil'*rho_stack))`.

### Local helper functions

- Line 80: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv and zeeman-liouv formalisms.')`.

## Parameters / inputs

- parameters.sweep nuclear frequency sweep width, Hz
- parameters.npoints number of fid points to be computed
- parameters.tau stimulated echo time, seconds
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -free induction decay whose Fourier transform is the
- Mims ENDOR signal

## Implementation structure

- Mims ENDOR pulse sequence with ideal hard pulses. Syntax:
- fid=endor_mims(spin_system,parameters,H,R,K)
- parameters.sweep nuclear frequency sweep width, Hz
- parameters.npoints number of fid points to be computed
- parameters.tau stimulated echo time, seconds
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- fid -free induction decay whose Fourier transform is the
- Mims ENDOR signal
- Move into adjoint representation if needed
- Consistency check

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `sim2liouv()`, `grumble()`, `state()`, `operator()`, `step()`, `evolution()`, `coherence()`, `transpose()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`.
