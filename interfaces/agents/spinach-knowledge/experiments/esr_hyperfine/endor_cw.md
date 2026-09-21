# experiments/esr_hyperfine/endor_cw.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/esr_hyperfine/endor_cw.m`
- Signature: `fid=endor_cw(spin_system,parameters,H,R,K)`
- Total lines: 97

## Purpose

Fast approximate simulation of isotropic continuous-wave ENDOR pulse sequence -essentially an NMR spectrum weighted by hyper- fine couplings is recorded. Syntax: fid=endor_cw(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Hyperfine ESR experiment implementations. These sequences probe coupled electron-nuclear dynamics through ENDOR or HYSCORE-type manipulations of coherence pathways.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters.sweep nuclear frequency sweep width, Hz
- parameters.npoints number of FID points to be computed
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -free induction decay whose Fourier transform
- approximates a CW ENDOR spectrum

## Implementation structure

- Fast approximate simulation of isotropic continuous-wave ENDOR
- pulse sequence -essentially an NMR spectrum weighted by hyper-
- fine couplings is recorded. Syntax:
- fid=endor_cw(spin_system,parameters,H,R,K)
- parameters.sweep nuclear frequency sweep width, Hz
- parameters.npoints number of FID points to be computed
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- fid -free induction decay whose Fourier transform
- approximates a CW ENDOR spectrum
- Move into adjoint representation if needed

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `sim2liouv()`, `grumble()`, `operator()`, `cellfun()`, `strncmp()`, `state()`, `step()`, `evolution()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`.
