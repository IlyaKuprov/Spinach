# kernel/utilities/adelim.m

- Signature: `[L,R]=adelim(spin_system,L,fast_idx,slow_idx)`

## Purpose

Perform adiabatic elimination in Liouville space, following Section 6.1 of Kuprov's book.

## Mathematical content

The input basis is partitioned into fast and slow indices. The function forms the four Liouvillian blocks `L00`, `L01`, `L10`, and `L11`; it returns the slow block `L=L00` and the additional relaxation superoperator `R=1i*L01*(L11\L10)`. The fast subsystem is expected to be dissipative.

## Parameters / inputs

- `spin_system` - Spinach spin-system structure; the supported formalisms are `sphten-liouv` and `zeeman-liouv`.
- `L` - square Liouvillian in the selected Liouville-space formalism.
- `fast_idx` - indices of basis states involving the fast subsystem.
- `slow_idx` - indices involving only the slow subsystem. In `sphten-liouv`, basis states can be attributed to individual spins; in `zeeman-liouv`, the caller must provide index sets meaningful in the Zeeman basis.

The fast and slow index sets must be disjoint and together cover the Liouvillian basis.

## Outputs

- `L` - projection of the original Liouvillian into the slow subspace, retaining its existing coherent and dissipative dynamics.
- `R` - additional relaxation superoperator from eliminating the fast subspace.
