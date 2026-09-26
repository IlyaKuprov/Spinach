# experiments/acquire.m

- Signature: `fid=acquire(spin_system,parameters,H,R,K)`

## Purpose

Simple forward time evolution with signal acquisition. Syntax: fid=acquire(spin_system,parameters,H,R,K)

## Physical / mathematical content

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Parameters / inputs

- parameters.sweep sweep width, Hz
- parameters.npoints number of points in the FID
- parameters.rho0 initial state
- parameters.coil detection state
- parameters.decouple spins to decouple, e.g. {'15N','13C'}
- parameters.homodec_oper operator to add to the Liouvillian at
- the detection stage
- parameters.homodec_pwr power coefficient for the operator, Hz
- parameters.dead_time the system will be evolved for this
- time (seconds) before the signal
- acquisition begins
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -free induction decay as seen by the state specified
- in parameters parameters.coil

## Implementation structure

- Simple forward time evolution with signal acquisition. Syntax:
- fid=acquire(spin_system,parameters,H,R,K)
- parameters.sweep sweep width, Hz
- parameters.npoints number of points in the FID
- parameters.rho0 initial state
- parameters.coil detection state
- parameters.decouple spins to decouple, e.g. {'15N','13C'}
- parameters.homodec_oper operator to add to the Liouvillian at
- the detection stage
- parameters.homodec_pwr power coefficient for the operator, Hz
- parameters.dead_time the system will be evolved for this
- time (seconds) before the signal
