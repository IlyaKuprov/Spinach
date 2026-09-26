# experiments/spin_chem/rydmr.m

- Signature: `A=rydmr(spin_system,parameters,H,R,K)`

## Purpose

Singlet-singlet RYDMR experiment using the full kinetics superoper- ator -computes the singlet yield of a radical pair recombination reaction. Syntax: A=rydmr(spin_system,parameters,H,R,K) where H is the Hamiltonian commutation superoperator in zero ex- ternal field, R is the relaxation superoperator and K is the che- mical kinetics superoperator. Parameters: parameters.tol -BICG solver tolerance, 1e-2 is generally g

## Physical / mathematical content

- Spin-chemistry experiment implementations. These routines couple spin evolution to chemical kinetics, radical-pair recombination, exchange, and spin-selective reaction channels.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Outputs

- A -fractional singlet yield

## Implementation structure

- Singlet-singlet RYDMR experiment using the full kinetics superoper-
- ator -computes the singlet yield of a radical pair recombination
- reaction. Syntax:
- A=rydmr(spin_system,parameters,H,R,K)
- where H is the Hamiltonian commutation superoperator in zero ex-
- ternal field, R is the relaxation superoperator and K is the che-
- mical kinetics superoperator. Parameters:
- parameters.tol - BICG solver tolerance,
- 1e-2 is generally good
- A -fractional singlet yield
- Check consistency
- Get the two-electron singlet state
