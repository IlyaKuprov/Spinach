# kernel/utilities/sinkhole.m

- Signature: `L=sinkhole(spin_system,L,states)`

## Purpose

Turns the specified states into sinkholes --any population reaching them will be summed up and stored forever in a frozen state. This is useful for state space restriction diagnostics. Syntax: L=sinkhole(spin_system,L,states)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- L -Liovillian matrix
- states -a vector of integers specifying the
- numbers of the states to be set up
- as sinkholes
- Output:
- L -updated Liouvillian matrix
- Note: this functionality is only available in sphten-liouv formalism.

## Implementation structure

- Turns the specified states into sinkholes --any population reaching
- them will be summed up and stored forever in a frozen state. This is
- useful for state space restriction diagnostics. Syntax:
- L=sinkhole(spin_system,L,states)
- L -Liovillian matrix
- states -a vector of integers specifying the
- numbers of the states to be set up
- as sinkholes
- Output:
- L -updated Liouvillian matrix
- Note: this functionality is only available in sphten-liouv formalism.
- Check consistency
