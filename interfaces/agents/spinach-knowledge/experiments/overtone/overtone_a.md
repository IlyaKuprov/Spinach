# experiments/overtone/overtone_a.m

- Signature: `spectrum=overtone_a(spin_system,parameters,H,R,K)`

## Purpose

Overtone signal acquisition experiment in the frequency domain. Syntax: spectrum=overtone_a(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Overtone experiment implementations. These routines excite or detect high-order quadrupolar transitions and therefore combine non-secular quadrupolar terms, MAS or field effects, and specialised detection pathways.

## Numerical / algorithmic content

## Parameters / inputs

- parameters.spins overtone-active nucleus, specified as a
- single-element cell array
- parameters.sweep vector with two elements giving
- the spectrum frequency extents
- in Hz around the overtone frequency
- parameters.npoints number of points in the spectrum
- parameters.rho0 initial state
- parameters.coil detection state
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- spectrum -the spectrum of the system with the specified
- starting state detected on the specified coil
- state within the frequency interval requested
- Note: relaxation must be present in the system dynamics, or the matrix
- inversion operation in the slowpass call would fail. The relaxa-
- tion superoperator R must *not* be thermalised.

## Implementation structure

- Overtone signal acquisition experiment in the frequency domain. Syntax:
- spectrum=overtone_a(spin_system,parameters,H,R,K)
- parameters.spins overtone-active nucleus, specified as a
- single-element cell array
- parameters.sweep vector with two elements giving
- the spectrum frequency extents
- in Hz around the overtone frequency
- parameters.npoints number of points in the spectrum
- parameters.rho0 initial state
- parameters.coil detection state
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
