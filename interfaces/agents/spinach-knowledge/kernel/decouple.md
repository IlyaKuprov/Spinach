# kernel/decouple.m

- Signature: `[L,rho]=decouple(spin_system,L,rho,spins)`

## Purpose

Removes interactions and state components involving the specified spins. Those spins do not contribute to the dynamics until the Liouvillian is rebuilt. The operation is the analytical equivalent of a perfect decoupling pulse sequence on the selected spins.

## Physical / mathematical content

- Requires `sphten-liouv`, `zeeman-liouv`, or `zeeman-hilb` formalism.
- In `sphten-liouv`, basis states involving a selected spin are zeroed.
- In the Zeeman formalisms, the selected spins are projected onto their identity components, because spin involvement is not diagonal in the Zeeman basis.
- In `zeeman-hilb`, the Hamiltonian and density matrices are stretched into Liouville space, projected, and folded back; each selected-spin Hamiltonian factor is replaced by its identity-component average.
- The Liouville-space formalisms support Fokker–Planck direct products; the operation is extended across the spatial subproblem.

## Parameters / inputs

- `L`: Liouvillian superoperator, or the Hamiltonian in `zeeman-hilb`; may be empty.
- `rho`: state vector or horizontal stack, or a density matrix/stack in `zeeman-hilb`; may be empty.
- `spins`: selected spins specified by names (for example, `{'13C','1H'}`) or indices (for example, `[1 2]`).

## Outputs

- `rho`: state(s) with components involving the selected spins removed.
- `L`: Liouvillian with terms involving the selected spins removed.
