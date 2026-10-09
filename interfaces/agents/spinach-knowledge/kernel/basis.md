# kernel/basis.m

Canonical implementation: [kernel/basis.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/basis.m)

- Signature: `spin_system=basis(spin_system,bas)`

## Purpose and representation

This is the mandatory basis-preparation step after [create.m](create.md). It stores `bas` in `spin_system.bas`, builds the basis descriptor for the requested formalism, reports basis information, and calls `symmetry`.

The accepted formalism strings are `zeeman-hilb`, `zeeman-wavef`, `zeeman-liouv`, and `sphten-liouv`. Every input field except `formalism` is a cell with exactly one entry per substance, including single-substance systems; no scalar broadcast is performed. `basis`, not the two-argument `create`, validates this contract.

`bas.nsubst` counts substances. `bas.basis{n}` is a sparse spherical-tensor descriptor with columns in the order of the global spin list `chem.parts{n}`; its first row is the unit. For other formalisms this cell is empty. `nstates(n)` is the descriptor row count for spherical tensors, the local Hilbert dimension for Hilbert/wavefunction formalisms, and its square for Zeeman Liouville space. A spin-free substance has one state. `offsets=[0;cumsum(nstates)]` addresses the direct sum.

Projection and correlation summaries are `tot_proj{n}` and `tot_cord{n}`. `sym_fact(n)` holds `irr_dimensions` and `irr_projectors`; projector columns are symmetry-adapted vectors in the local descriptor coordinates, as in the existing SALC algorithm. The local SALC service itself uses a one-cell descriptor and returns `sym_fact`, without an intermediate legacy `irrep` structure. No global descriptor or `bas.irrep` is returned. Single-substance Zeeman symmetry uses temporary stock ket/bra index labels for the SALC calculation and stores the resulting projectors in `sym_fact(1)`; the compiled Zeeman descriptor cell remains empty. Multi-substance Zeeman symmetry raises `Spinach:basis:segmentedZeeman` pending its implementation.

## Approximation and graph options

`bas.approximation` is required. The Zeeman Hilbert and Zeeman Liouville formalisms require `none` in each cell; `sphten-liouv` accepts `none`, `IK-0`, `IK-1`, `IK-2`, `IK-DNP`, or `IK-SBS`.

- `IK-1` and `IK-2` require `bas.prox_level`; the positive integer depth cannot exceed the number of spins. `IK-0`, `IK-1`, `IK-DNP`, and `IK-SBS` require `bas.inter_level`. For `IK-0`/`IK-1` it is a positive integer no larger than the spin count. For `IK-DNP` it is a three-entry positive-integer vector bounded respectively by electron count, spin count, and nucleus count. For `IK-SBS` its three entries are bounded respectively by bosonic-mode count, total particle count, and spin count.
- `IK-1`, `IK-2`, and `IK-SBS` require `bas.connectivity`, either `scalar_couplings` or `full_tensors`. The scalar option uses `abs(trace(tensor)/3)`; the full-tensor option uses the matrix 2-norm `norm(tensor,2)`. Couplings are included in the graph when their norm exceeds `2*pi*spin_system.tols.inter_cutoff`; the configured cutoff is reported in Hz. The graph is made reciprocal. `IK-DNP` constructs separate electron-electron, electron-nuclear, and nuclear-nuclear connectivity graphs; `IK-SBS` also builds mode connectivity from bosonic couplings and coupling derivatives.
- `IK-1` and `IK-2` are restricted to spin-only systems. `IK-DNP` requires both electrons and nuclei and rejects other component types. `IK-SBS` requires both spin and bosonic modes (types `C`, `V`, or `T`).
- In the spherical-tensor path, each substance in `spin_system.chem.parts` is handled using its corresponding interaction/connectivity block. Subgraphs and proximity neighborhoods determine candidate states; a unit state is retained, filters are applied, duplicate states are removed, and each local descriptor is sorted independently. `spin_system.bas.tot_proj{n}` and `spin_system.bas.tot_cord{n}` store each row's total projection quantum number and correlation order.

## State filters

The optional `bas.projections`, `bas.longitudinal`, and `bas.zero_quantum` settings are available only for `sphten-liouv`. Each is a cell array with one entry per chemical substance; `projections` entries are real integer row vectors (or empty), while `longitudinal` and `zero_quantum` entries are cells naming spins by isotope string or spin-index vector. `bas.projections` selects coherence orders while retaining the unit state. `bas.zero_quantum` selects zero total projection over its specified spins. `bas.longitudinal` keeps only longitudinal single-spin states on its specified spins. Optional `bas.manual{n}` is a logical table with one column per local spin; each included row specifies an additional subgraph. Symmetry spin indices are also local; numeric longitudinal and zero-quantum filter indices remain global. Empty entries for inapplicable depth/connectivity fields allow different approximations in different substances. `space_level` is an alternative spelling of `prox_level`; specifying both is rejected. These settings act on the spherical-tensor basis descriptor, not on the Hilbert Zeeman table.

## Cache-related state

For `sphten-liouv`, the routine prepares `spin_system.bas.lpst` and `spin_system.bas.rpst` cells indexed by spin multiplicity, loading or computing product tables through `ist_product_table` for each present multiplicity other than 1. It stores an MD5 hash of all descriptor blocks, `nstates`, and `chem.parts` in `spin_system.bas.basis_hash` for later cache use. This function does not itself describe persistent cache storage or cache invalidation.

## Parameters / output

- `spin_system` - primary structure produced by [create.m](create.md), including component multiplicities, types, isotope identities, interactions, and chemical-substance mapping.
- `bas` - basis specification; required fields depend on formalism and approximation as summarised above.
- Output `spin_system` - updated with `bas`, the selected basis and related metadata, and symmetry treatment.

## Reference

The source points to [10.1063/1.3624564](http://link.aip.org/link/doi/10.1063/1.3624564) for discussion of basis-set selection.

[Spin Dynamics Wiki: basis.m](https://spindynamics.org/wiki/index.php?title=basis.m)

Each substance starts with fresh local options. In particular, deriving `prox_level` from a nonempty `space_level` entry cannot leak that derived value into a later empty entry.

Before compilation, the union of `chem.parts` must equal every global spin index; an omitted spin raises `Spinach:basis:incompletePartition`. Empty spin-free substances remain valid when the other parts cover all spins.

Wavefunction direct-sum storage uses the local Hilbert dimensions, including scalar spin-free blocks. Nonempty reaction records raise `Spinach:basis:wavefunctionChemistry`: chemical reactions are not supported in zeeman-wavef formalism.

Legacy global `bas.basis` matrices and `bas.irrep` fields are rejected at this entry point with named errors pointing to per-substance `bas.basis{n}`/`bas.offsets` and `bas.sym_fact(n)` symmetry data. Compiled structures remain ordinary MATLAB structs; arbitrary external dot reads are not intercepted.
