# kernel/basis.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/basis.m`
- Signature: `spin_system=basis(spin_system,bas)`
- Total lines: 843

## Purpose

Basis set control. This is the second mandatory function (after create.m) that must be called in every calculation to build spin_system data structure. Syntax: spin_system=basis(spin_system,bas)

## Physical / mathematical content

- In `sphten-liouv` formalism, subgraphs are generated separately for each chemical substance listed in `spin_system.chem.parts`, using the connectivity and proximity information of that substance only; the correlation levels `bas.inter_level` and `bas.prox_level` are clipped to the spin count of each substance.
- The state filters `bas.projections`, `bas.longitudinal`, and `bas.zero_quantum` are cell arrays with one element per substance; an empty element means no filter for that substance. The zero-quantum filter keeps states whose total projection over the union of the listed spins is zero.
- The per-substance state lists are merged into one global basis with a single unit state and sorted lexicographically; `spin_system.bas.tot_proj` and `spin_system.bas.tot_cord` hold the total projection quantum number and the correlation order of each basis state.

## Numerical / algorithmic content

- Subgraph generation uses `dfpt` on the substance blocks of the connectivity and proximity matrices; empty, identical, and enclosed subgraphs are removed with `unique` and `prune_subgraphs` before the descriptor is built.
- The descriptor of each subgraph is built densely in the direct product order with `repelem`/`repmat`, filtered, and embedded into the full spin index as a sparse array; duplicate states across subgraphs are removed with `unihash`, and the basis is sorted with `sortrows`, distributed for large bases.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 46-47: Show the banner; implemented by `banner(spin_system,'basis_banner')`.
- Lines 49-50: Check the input; implemented by `grumble(spin_system,bas)`.
- Lines 52-53: Store the settings; implemented by `spin_system.bas=bas`.
- Lines 55-57: Find electrons and nuclei; implemented by `e_idx=cellfun(@iselectron,spin_system.comp.isotopes)`.
- Lines 59-60: Report back to the user; implemented by `summary_basis_opts(spin_system)`.
- Lines 62-63: Process spherical tensor basis sets; implemented by `if strcmp(spin_system.bas.formalism,'sphten-liouv')`.
- Lines 65-68: Disallow spherical tensor basis sets for large multiplicities; implemented by `if any(spin_system.comp.mults>16,'all')`.
- Lines 70-71: Count chemical substances; implemented by `nsubst=numel(spin_system.chem.parts)`.
- Lines 73-111: Run connectivity analysis for IK-DNP basis set; builds `ee_conmatrix`, `en_conmatrix`, and `nn_conmatrix` from the coupling tensors above `spin_system.tols.inter_cutoff`.
- Lines 113-168: Run connectivity analysis for IK-1,2 basis sets; builds `spin_system.inter.conmatrix` from `bas.connectivity` and the bosonic mode channels, symmetrises it and `spin_system.inter.proxmatrix`, and reports their densities.
- Lines 170-174: Build state lists for individual spins; implemented by `spin_state_lists{n}=(0:(spin_system.comp.mults(n)^2-1))'`.
- Lines 176-203: Apply longitudinal filters; for each substance, resolves `bas.longitudinal{s}` entries to spins of that substance and removes states with `M~=0` from their lists.
- Lines 205-206: Compute subspace dimensions for individual spins; implemented by `spin_dims=cellfun(@numel,spin_state_lists)`.
- Lines 211-348: Loop over chemical substances; generates coupling and proximity subgraphs in the spin index of the substance (`none`, `IK-0`, `IK-1`, `IK-2`, `IK-DNP`), adds the rows of `bas.manual` that belong to the substance, removes spin-zero particles, empty, identical, and enclosed subgraphs, and embeds the result into the full spin index.
- Lines 350-352: Merge the subgraph lists of all substances; implemented by `subgraphs=vertcat(subgraphs{:})`.
- Lines 354-382: Resolve zero-quantum filter spins for each substance; implemented by `zq_spins(s,spins_in_question)=true()`.
- Lines 384-386: Balance the subgraph list; implemented by `randperm`.
- Lines 388-459: Populate the basis descriptor array; `parfor` over subgraphs builds the direct product descriptor, applies the coherence order, zero-quantum, and IK-DNP inter-nuclear filters, and embeds it into the full spin index with `sparse`.
- Lines 461-463: Eliminate redundant states using a hash table; implemented by `basis_spec=unihash(basis_spec)`.
- Lines 465-478: Sort the basis explicitly; implemented by `sortrows`, distributed with `distrib_dim` for large bases.
- Lines 484-488: Total projection quantum number and correlation order of each state; implemented by `[L,M]=lin2lm(spin_system.bas.basis)`.
- Lines 490-494: Report on chemical species; state counts per substance.
- Lines 496-500: Print the summary and run the symmetry treatment; implemented by `summary_basis(spin_system)` and `symmetry(spin_system,bas)`.
- Lines 504-529: Process Hilbert space Zeeman basis; Zeeman index table for `zeeman-hilb` and `zeeman-wavef`.
- Lines 531-558: Process Liouville space Zeeman basis; ket and bra index tables for `zeeman-liouv`.
- Lines 560-587: Preload Lie algebra structure tables; implemented by `ist_product_table`.
- Lines 589-593: Hash the basis descriptor for caching tools later; implemented by `md5_hash(spin_system.bas.basis)`.

### Control flow inferred from the code

- Line 63: conditional branch on `strcmp(spin_system.bas.formalism,'sphten-liouv')`.
- Line 74: conditional branch on `strcmp(spin_system.bas.approximation,'IK-DNP')`.
- Line 114: conditional branch on `ismember(spin_system.bas.approximation,{'IK-1','IK-2'})`.
- Line 177: conditional branch on `isfield(bas,'longitudinal')`.
- Line 212: `for` loop over `s=1:nsubst`.
- Line 220: dispatches on `spin_system.bas.approximation`; cases `'none'`, `'IK-0'`, `'IK-1'`, `'IK-2'`, `'IK-DNP'`.
- Line 310: conditional branch on `isfield(bas,'manual')`; a row of `bas.manual` that spans two substances is an error.
- Line 355: conditional branch on `isfield(bas,'zero_quantum')`.
- Line 391: `parfor` loop over `n=1:size(subgraphs,1)`.
- Line 505: conditional branch on `ismember(spin_system.bas.formalism,{'zeeman-hilb','zeeman-wavef'})`.
- Line 532: conditional branch on `strcmp(spin_system.bas.formalism,'zeeman-liouv')`.

### Local helper functions

- Line 598: `grumble()` — `function grumble(spin_system,bas)`. Checks `bas.formalism`, `bas.approximation`, `bas.connectivity`, `bas.inter_level`, `bas.prox_level`, `bas.manual`, `bas.projections`, `bas.longitudinal`, and `bas.zero_quantum`; the three filters must be cell arrays with one element per chemical substance and are only available for `sphten-liouv`.

## Parameters / inputs

- spin_system - primary Spinach data structure, the output of create.m function
- bas - basis set specification structure described in detail in the online manual

## Outputs

- spin_system - primary Spinach data structure, updated with the basis set and related information
- Note: it is important to understand the factors that influence basis set selection in spin dynamics simulations - see our paper for further information on this subject.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `banner()`, `grumble()`, `cellfun()`, `summary_basis_opts()`, `report()`, `dfpt()`, `nchoosek()`, `prune_subgraphs()`, `unique()`, `lin2lm()`, `repelem()`, `sparse()`, `unihash()`, `distrib_dim()`, `sortrows()`, `summary_basis()`, `symmetry()`, `ist_product_table()`, `md5_hash()`.
