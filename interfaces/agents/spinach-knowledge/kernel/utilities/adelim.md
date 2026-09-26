# kernel/utilities/adelim.m

- Signature: `[L,R]=adelim(spin_system,L,fast_idx,slow_idx)`

## Purpose

Adiabatic elimination in Liouville space, implements Section 6.1 of Kuprov's book. Syntax: [L,R]=adelim(spin_system,L,fast_idx,slow_idx)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- L -Liouvillian in a Liouville space formalism,
- fast subbsystem must be dissipative
- fast_idx -a vector of integers specifying which
- states in the basis involve the fast
- subsystem in any way
- slow_idx -a vector of integers specifying which
- states in the basis only involve the
- slow subsystem

## Outputs

- L -projection of the original Liouvillian
- into the slow subspace, inheriting any
- coherent and dissipative dynamics that
- the user previously had there
- R -the extra relaxation superoperator on-
- ce the fast subspace is adiabatically
- eliminated
- Note: in sphten-liouv the basis states are attributable to
- individual spins; in zeeman-liouv the caller must
- supply index sets that are meaningful in the Zeeman
- basis of Liouville space.

## Implementation structure

- Adiabatic elimination in Liouville space, implements
- Section 6.1 of Kuprov's book. Syntax:
- [L,R]=adelim(spin_system,L,fast_idx,slow_idx)
- L -Liouvillian in a Liouville space formalism,
- fast subbsystem must be dissipative
- fast_idx -a vector of integers specifying which
- states in the basis involve the fast
- subsystem in any way
- slow_idx -a vector of integers specifying which
- states in the basis only involve the
- slow subsystem
- L -projection of the original Liouvillian
