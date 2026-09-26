# experiments/esr_hyperfine/endor_mims.m

- Signature: `fid=endor_mims(spin_system,parameters,H,R,K)`

## Purpose

Mims ENDOR pulse sequence with ideal hard pulses. Syntax: fid=endor_mims(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Hyperfine ESR experiment implementations. These sequences probe coupled electron-nuclear dynamics through ENDOR or HYSCORE-type manipulations of coherence pathways.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

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
