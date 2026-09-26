# experiments/traject.m

- Signature: `traj=traject(spin_system,parameters,H,R,K)`

## Purpose

Simple forward time evolution trajectory. Syntax: traj=traject(spin_system,parameters,H,R,K)

## Physical / mathematical content

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Parameters / inputs

- parameters.sweep sweep width, Hz
- parameters.npoints number of points in the trajectory
- parameters.rho0 initial state
- parameters.decouple spins to decouple, e.g. {'15N','13C'}
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- traj -system trajectory, a bookshelf stack of state vectors

## Implementation structure

- Simple forward time evolution trajectory. Syntax:
- traj=traject(spin_system,parameters,H,R,K)
- parameters.sweep sweep width, Hz
- parameters.npoints number of points in the trajectory
- parameters.rho0 initial state
- parameters.decouple spins to decouple, e.g. {'15N','13C'}
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- traj -system trajectory, a bookshelf stack of state vectors
- Check consistency
- Compose Liouvillian
