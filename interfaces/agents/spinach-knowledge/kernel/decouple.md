# kernel/decouple.m

- Signature: `[L,rho]=decouple(spin_system,L,rho,spins)`

## Purpose

Obliterates all interactions and populations in the subspace of states that involve the specified spins in any way. The specified spins would not contribute to the system dynamics until the Liouvillian is rebuilt from scratch. Syntax: [L,rho]=decouple(spin_system,L,rho,spins)

## Physical / mathematical content

- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

## Parameters / inputs

- L -Liouvillian superoperator or, in zeeman-hilb,
- the Hamiltonian; this may be left empty
- rho -state vector or a horizontal stack thereof or,
- in zeeman-hilb, a density matrix or a horizon-
- tal stack thereof; this may be left empty
- spins -spins to be wiped, specified either by name, e.g.
- {'13C','1H'}, or by a list of numbers, e.g. [1 2]

## Outputs

- rho -state vector(s) with all populations of the
- states involving the target spins set to zero
- L -Liouvillian superoperator with all terms in-
- volving the target spins set to zero
- Note: this function is an analytical equivalent of a perfect decoup-
- ling pulse sequence on the specified spins.
- Note: this function requires sphten-liouv, zeeman-liouv, or zeeman-
- hilb formalism; Fokker-Planck direct products are supported in
- the Liouville space formalisms. In the Zeeman formalisms the
- operation is an exact projection onto the subspace where the
- decoupled spins carry only their identity component, because
- spin involvement is not diagonal in the Zeeman basis. In
- zeeman-hilb, the Hamiltonian and the density matrices are
- stretched into Liouville space, projected there, and folded
- back; this replaces every decoupled-spin factor of the Hamil-
- tonian by its identity component average.

## Implementation structure

- Obliterates all interactions and populations in the subspace of states
- that involve the specified spins in any way. The specified spins would
- not contribute to the system dynamics until the Liouvillian is rebuilt
- from scratch. Syntax:
- [L,rho]=decouple(spin_system,L,rho,spins)
- L -Liouvillian superoperator or, in zeeman-hilb,
- the Hamiltonian; this may be left empty
- rho -state vector or a horizontal stack thereof or,
- in zeeman-hilb, a density matrix or a horizon-
- tal stack thereof; this may be left empty
- spins -spins to be wiped, specified either by name, e.g.
- {'13C','1H'}, or by a list of numbers, e.g. [1 2]
