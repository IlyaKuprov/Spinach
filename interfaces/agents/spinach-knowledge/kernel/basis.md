# kernel/basis.m

Canonical implementation: [kernel/basis.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/basis.m)

- Signature: `spin_system=basis(spin_system,bas)`

## Purpose and representation

This is the mandatory basis-preparation step after [create.m](create.md). It stores `bas` in `spin_system.bas`, builds the basis descriptor for the requested formalism, reports basis information, and calls `symmetry`.

The accepted formalism strings are `zeeman-hilb`, `zeeman-wavef`, `zeeman-liouv`, and `sphten-liouv`. Let `d=prod(spin_system.comp.mults)` be the Hilbert-space dimension and `N=spin_system.comp.nspins`.

- For `zeeman-hilb` and `zeeman-wavef`, the Zeeman-label table has `d` rows and `N` spin columns.
- For `zeeman-liouv`, the table has `d^2` rows and `2*N` columns: the implementation combines Hilbert Zeeman labels into ket/bra index tables. It reports dimension `d^2` for superoperators and state vectors.
- For `sphten-liouv`, rows are the generated, filtered spherical-tensor states; the table is merged, sorted, and accompanied by per-state total projection and correlation-order arrays. The state count depends on approximation and filters, rather than being fixed by `d`. Multiplicities above 16 are rejected for this formalism.

## Approximation and graph options

`bas.approximation` is required. The Zeeman Hilbert and Zeeman Liouville formalisms require `none`; `sphten-liouv` accepts `none`, `IK-0`, `IK-1`, `IK-2`, `IK-DNP`, or `IK-SBS`.

- `IK-1` and `IK-2` require `bas.prox_level`; the positive integer depth cannot exceed the number of spins. `IK-0`, `IK-1`, `IK-DNP`, and `IK-SBS` require `bas.inter_level`. For `IK-0`/`IK-1` it is a positive integer no larger than the spin count. For `IK-DNP` it is a three-entry positive-integer vector bounded respectively by electron count, spin count, and nucleus count. For `IK-SBS` its three entries are bounded respectively by bosonic-mode count, total particle count, and spin count.
- `IK-1`, `IK-2`, and `IK-SBS` require `bas.connectivity`, either `scalar_couplings` or `full_tensors`. The scalar option uses `abs(trace(tensor)/3)`; the full-tensor option uses the matrix 2-norm `norm(tensor,2)`. Couplings are included in the graph when their norm exceeds `2*pi*spin_system.tols.inter_cutoff`; the configured cutoff is reported in Hz. The graph is made reciprocal. `IK-DNP` constructs separate electron-electron, electron-nuclear, and nuclear-nuclear connectivity graphs; `IK-SBS` also builds mode connectivity from bosonic couplings and coupling derivatives.
- `IK-1` and `IK-2` are restricted to spin-only systems. `IK-DNP` requires both electrons and nuclei and rejects other component types. `IK-SBS` requires both spin and bosonic modes (types `C`, `V`, or `T`).
- In the spherical-tensor path, each substance in `spin_system.chem.parts` is handled using its corresponding interaction/connectivity block. Subgraphs and proximity neighborhoods determine candidate states; a unit state is retained, filters are applied, duplicate states are removed, and the global descriptor is sorted. `spin_system.bas.tot_proj` and `spin_system.bas.tot_cord` store each row's total projection quantum number and correlation order.

## State filters

The optional `bas.projections`, `bas.longitudinal`, and `bas.zero_quantum` settings are available only for `sphten-liouv`. Each is a cell array with one entry per chemical substance; `projections` entries are real integer row vectors (or empty), while `longitudinal` and `zero_quantum` entries are cells naming spins by isotope string or spin-index vector. `bas.projections` selects coherence orders while retaining the unit state. `bas.zero_quantum` selects zero total projection over its specified spins. `bas.longitudinal` keeps only longitudinal single-spin states on its specified spins. Optional `bas.manual` is a logical or numeric table with one column per spin; each included row specifies an additional subgraph. These settings act on the spherical-tensor basis descriptor, not on the Hilbert Zeeman table.

## Cache-related state

For `sphten-liouv`, the routine prepares `spin_system.bas.lpst` and `spin_system.bas.rpst` cells indexed by spin multiplicity, loading or computing product tables through `ist_product_table` for each present multiplicity other than 1. If `op_cache` or `ham_cache` is enabled in `spin_system.sys.enable`, it stores `md5_hash(spin_system.bas.basis)` in `spin_system.bas.basis_hash` for later cache use. This function does not itself describe persistent cache storage or cache invalidation.

## Parameters / output

- `spin_system` - primary structure produced by [create.m](create.md), including component multiplicities, types, isotope identities, interactions, and chemical-substance mapping.
- `bas` - basis specification; required fields depend on formalism and approximation as summarised above.
- Output `spin_system` - updated with `bas`, the selected basis and related metadata, and symmetry treatment.

## Reference

The source points to [10.1063/1.3624564](http://link.aip.org/link/doi/10.1063/1.3624564) for discussion of basis-set selection.

[Spin Dynamics Wiki: basis.m](https://spindynamics.org/wiki/index.php?title=basis.m)
