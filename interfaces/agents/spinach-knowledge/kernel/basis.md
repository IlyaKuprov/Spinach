# kernel/basis.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/basis.m`
- Signature: `spin_system=basis(spin_system,bas)`
- Total lines: 1018

## Purpose

Basis set control. This is the second mandatory function (after create.m) that must be called in every calculation to build spin_system data structure. Syntax: spin_system=basis(spin_system,bas)

## Physical / mathematical content

- In `sphten-liouv` formalism, subgraphs are generated separately for each chemical substance listed in `spin_system.chem.parts`, using the connectivity and proximity information of that substance only; the correlation levels `bas.inter_level` and `bas.prox_level` are clipped to the spin count of each substance.
- The state filters `bas.projections`, `bas.longitudinal`, and `bas.zero_quantum` are cell arrays with one element per substance; an empty element means no filter for that substance. The zero-quantum filter keeps states whose total projection over the union of the listed spins is zero.
- `IK-1` and `IK-2` are spin-only: the grumbler refuses systems with bosonic modes, and their coupling graph is built from the spin-spin coupling tensors alone. `IK-SBS` (spin-boson systems) builds its own graph from the spin-spin tensors under the same `bas.connectivity` norm plus all bosonic mode couplings above `tols.inter_cutoff` (pairwise `exchange`, `dispersive`, `kerr`, `longitudinal` channels by 2-norm; `coupling_mod` and `zeeman_mod` entries link the modes of the outer cell to the modulated spin pairs and spins of the inner block), splits it into boson-boson, spin-boson, and spin-spin graphs (bosonic modes are the `C`, `V`, and `T` particles of `spin_system.comp.types`), traces each with `dfpt` to its own level in `bas.inter_level=[bb sb ss]`, and merges the three subgraph lists; inside the resulting subgraphs, pure boson-boson correlations above the first level and pure spin-spin correlations above the third level are dropped, mirroring the inter-nuclear filter of `IK-DNP`. Both spins and modes must be present, and `bas.connectivity` is required.
- The per-substance state lists are merged into one global basis with a single unit state and sorted lexicographically; `spin_system.bas.tot_proj` and `spin_system.bas.tot_cord` hold the total projection quantum number and the correlation order of each basis state.

## Numerical / algorithmic content

- Subgraph generation uses `dfpt` on the substance blocks of the connectivity and proximity matrices; empty, identical, and enclosed subgraphs are removed with `unique` and `prune_subgraphs` before the descriptor is built.
- The descriptor of each subgraph is built densely in the direct product order with `repelem`/`repmat` in the smallest signed integer class that holds every single-spin state index of the system (`min_int_type` of the largest multiplicity squared minus one: `int8` up to multiplicity 11, `int16` above), filtered through the integer-aware `lin2lm`, and embedded into the full spin index as a sparse double array.
- Duplicate states across subgraphs are removed and the basis is sorted lexicographically in one `unique(...,'rows')` call; when a dense integer copy of the merged descriptor takes fewer bytes than its sparse form (one byte per element against sixteen per non-zero), the call runs on that dense copy and the result is converted back, otherwise it runs on the sparse matrix. The stored `spin_system.bas.basis` is sparse double in either case, so every consumer sees the same descriptor as before.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.

## Parameters / inputs

- spin_system - primary Spinach data structure, the output of create.m function
- bas - basis set specification structure described in detail in the online manual

## Outputs

- spin_system - primary Spinach data structure, updated with the basis set and related information
- Note: it is important to understand the factors that influence basis set selection in spin dynamics simulations - see our paper for further information on this subject.
