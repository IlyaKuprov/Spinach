# kernel/utilities/adelim.m

## Purpose

`adelim.m` performs adiabatic elimination in Liouville space, implementing Section 6.1 of Kuprov's book. It projects a Liouvillian into the slow subspace and returns the extra relaxation superoperator that appears once the fast subspace is adiabatically eliminated.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/adelim.m>

## Behaviour

The function is called as `[L,R]=adelim(spin_system,L,fast_idx,slow_idx)`.

1. A consistency check (`grumble`) validates the inputs: the formalism must be `sphten-liouv` or `zeeman-liouv`, `L` must be a square numeric matrix, `fast_idx` and `slow_idx` must be numeric vectors with no common elements, and the total number of indices must equal the dimension of `L`.
2. The unit operator `U` is obtained via `unit_oper(spin_system)`.
3. Projectors are formed as `P_slow=U(slow_idx,:)` and `P_fast=U(fast_idx,:)` (noted in the source as faster than indexing).
4. The Liouvillian is partitioned into blocks: `L01=P_slow*L*P_fast'`, `L10=P_fast*L*P_slow'`, `L11=P_fast*L*P_fast'`, and `L00=P_slow*L*P_slow'`.
5. Following Eq. 6.2 in Kuprov's book, the relaxation superoperator is computed as `R=1i*L01*(L11\L10)` and the projected Liouvillian is `L=L00`.

The header notes that in `sphten-liouv` the basis states are attributable to individual spins, while in `zeeman-liouv` the caller must supply index sets that are meaningful in the Zeeman basis of Liouville space. The fast subsystem must be dissipative.

## Inputs and outputs

Inputs:

- `spin_system` — spin system object.
- `L` — Liouvillian in a Liouville space formalism; the fast subsystem must be dissipative.
- `fast_idx` — vector of integers specifying which states in the basis involve the fast subsystem in any way.
- `slow_idx` — vector of integers specifying which states in the basis only involve the slow subsystem.

Outputs:

- `L` — projection of the original Liouvillian into the slow subspace, inheriting any coherent and dissipative dynamics that the user previously had there.
- `R` — the extra relaxation superoperator once the fast subspace is adiabatically eliminated.

## References

- Section 6.1 and Eq. 6.2 of Kuprov's book (as cited in the source header).
- Spinach Wiki page for `adelim.m`: <https://spindynamics.org/wiki/index.php?title=adelim.m>
- Source file: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/adelim.m>
