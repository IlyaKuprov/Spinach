# experiments/spin_chem/rydmr.m

- Signature: `A=rydmr(spin_system,parameters,H,R,K)`

## Purpose

Computes the fractional singlet yield for radical-pair recombination in a singlet-singlet RYDMR experiment.

## Physical / mathematical content

- Constructs and normalizes a two-electron singlet state. The source forms the kinetics Liouvillian `L=H+1i*R+1i*K` from the supplied Hamiltonian, relaxation, and chemical-kinetics superoperators.
- Computes the yield from the singlet projection of the linear-system solution, weighted by the first radical-pair reaction rate.

## Numerical / algorithmic content

- Solves the source-defined Liouville-space linear system with BICG, using `parameters.tol` as the solver tolerance; the source says `1e-2` is generally a good value. This is a stationary yield calculation, not a time-domain propagation.
- The code can move the calculation to the GPU when enabled and gathers the result for output.

## Outputs

- A -fractional singlet yield

## Implementation structure

- Checks the dimensions of `H`, `R`, and `K` and the tolerance parameter.
- Builds the singlet state and Liouvillian, obtains the BICG solution, evaluates the fractional yield, and returns it.
