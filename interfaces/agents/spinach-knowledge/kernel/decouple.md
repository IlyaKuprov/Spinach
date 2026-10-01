# kernel/decouple.m

Source: [kernel/decouple.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/decouple.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=decouple.m)

- Signature: `[L,rho]=decouple(spin_system,L,rho,spins)`

## Purpose

Removes selected-spin involvement from the supplied operator and/or state data. It is an analytical decoupling operation, not a pulse-sequence simulation; after it is applied, those spins do not contribute to the dynamics until the Liouvillian is rebuilt.

## Spin-system context

`create(sys,inter)` constructs the `spin_system` object consumed by this routine. `decouple` reads its formalism, basis, isotope labels, spin multiplicities/count, and Liouvillian cleanup tolerance; it does not assemble or recreate the spin system. See the sibling [create.m page](create.md).

## Inputs and outputs

- `spin_system`: the created Spinach spin-system structure. The supported formalisms are `sphten-liouv`, `zeeman-liouv`, and `zeeman-hilb`; Fokker–Planck direct products are supported in the Liouville-space formalisms.
- `L`: a square Liouvillian, or a Hamiltonian in `zeeman-hilb`; it may be empty. When both nonempty `L` and `rho` are supplied, the code requires the column count of `L` to equal the row count of `rho`.
- `rho`: a state vector or a stack of state vectors in Liouville space; in `zeeman-hilb`, a density matrix or a stack of density matrices. It may be empty. The spin-space basis size is taken from `size(spin_system.bas.basis,1)`; the routine treats additional Fokker–Planck coordinates as a direct-product factor.
- `spins`: selected spins as a cell array of isotope names (for example, `{'13C','1H'}`) or numeric spin indices (for example, `[1 2]`). An empty selection returns immediately.
- Outputs `L` and `rho` are modified only when requested as outputs and their corresponding inputs are nonempty.

## Projection performed

In `sphten-liouv`, the zero mask flags basis rows for which the sum of the selected-spin basis entries is nonzero. Those state components are zeroed in `rho`, and the matching rows and columns of `L` are zeroed. With a Fokker–Planck direct product, this mask is repeated across the spatial/orientational subspace. For polyadic Liouvillians, the same zero mask is applied as a diagonal projector on both sides, preserving unopened implicit cores rather than requiring indexed assignment.

In the Zeeman formalisms, spin involvement is not diagonal in the Zeeman basis. The code instead constructs a projector onto the selected spins' identity components, using the spin multiplicities and matrix-unit factors. For `zeeman-liouv`, the projector is extended over any Fokker–Planck coordinates; it is applied to both sides of `L` and to `rho`. For `zeeman-hilb`, the Hamiltonian and density-matrix stack are reshaped into Liouville-space columns, projected, and reshaped back. Thus each selected-spin Hamiltonian/state factor is reduced to its identity-component average, as described in the source comments.

After a requested nonempty `L` is projected, the routine calls `clean_up` with `spin_system.tols.liouv_zero`.

## Guards

The source rejects unsupported formalisms, nonsquare `L`, incompatible nonempty `L`/`rho` dimensions, unknown isotope labels, and numeric spin indices that are non-real, below 1, nonintegral, or above the system spin count. The call may omit either data argument by passing it empty.
