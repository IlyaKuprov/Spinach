# experiments/respiration.m

- Signature: `fid=respiration(spin_system,parameters,H,R,K)`

## Purpose

RESPIRATION cross-polarisation method described in the paper from the Aarhus group (http://dx.doi.org/10.1021/jz3000905). Syntax: fid=respiration(spin_system,parameters,H,R,K)

## Physical / mathematical content

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Parameters / inputs

- parameters.sweep sweep width, Hz
- parameters.npoints number of points in the FID
- parameters.rho0 initial state
- parameters.coil detection state
- parameters.nloops number of RESPIRATION loops
- parameters.theta the angle of the ideal pulse
- at the end of each loop
- parameters.rate RESPIRATION pulse train rate, Hz
- parameters.spins working spins, e.g. {'1H','13C'}
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -free induction decay as seen by the state specified
- in parameters parameters.coil

## Implementation structure

- RESPIRATION cross-polarisation method described in the paper from
- the Aarhus group (http://dx.doi.org/10.1021/jz3000905). Syntax:
- fid=respiration(spin_system,parameters,H,R,K)
- parameters.sweep sweep width, Hz
- parameters.npoints number of points in the FID
- parameters.rho0 initial state
- parameters.coil detection state
- parameters.nloops number of RESPIRATION loops
- parameters.theta the angle of the ideal pulse
- at the end of each loop
- parameters.rate RESPIRATION pulse train rate, Hz
- parameters.spins working spins, e.g. {'1H','13C'}
