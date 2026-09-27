# experiments/overtone/overtone_a.m

- Signature: `spectrum=overtone_a(spin_system,parameters,H,R,K)`

## Purpose

Frequency-domain overtone acquisition: the routine converts the requested sweep offsets to absolute frequencies around the overtone reference and delegates acquisition to `slowpass`.
## Physical / mathematical content

The overtone reference frequency is `-2*spin(parameters.spins{1})*spin_system.inter.magnet/(2*pi)`. The sweep is specified relative to this frequency.
## Numerical / algorithmic content

The function validates its inputs, computes the overtone reference frequency, changes `parameters.sweep` to `ovt_frq-parameters.sweep`, and calls `slowpass` with the adjusted parameters and supplied dynamics matrices.
## Parameters / inputs

- parameters.spins overtone-active nucleus, specified as a single-element cell array
- parameters.sweep vector with two elements giving the spectrum frequency extents in Hz around the overtone frequency
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

A thin wrapper around `slowpass`, preceded by the local `grumble` input validator.