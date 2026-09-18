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
- The descriptor of each subgraph is built densely in the direct product order with `repelem`/`repmat`, filtered, and embedded into the full spin index as a sparse array; duplicate states across subgraphs are removed with `unihash`, and the basis is sorted with `sortrows`, distributed for large bases.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 32-33: Show the banner; implemented by `banner(spin_system,'basis_banner')`.
- Lines 35-36: Check the input; implemented by `grumble(spin_system,bas)`.
- Lines 38-39: Store the settings; implemented by `spin_system.bas=bas`.
- Lines 41-43: Find electrons and nuclei; implemented by `e_idx=cellfun(@iselectron,spin_system.comp.isotopes)`.
- Lines 45-46: Find bosonic modes; implemented by `b_idx=ismember(spin_system.comp.types,{'C','V','T'});`
- Lines 48-49: Report back to the user; implemented by `summary_basis_opts(spin_system)`.
- Lines 51-54: Remind the user about the amplitude cut-off; implemented by `report(spin_system,['coupling tensors with norm below ' ...`
- Lines 56-655: Process spherical tensor basis sets; implemented by `if strcmp(spin_system.bas.formalism,'sphten-liouv')`.
- Lines 59-62: Disallow spherical tensor basis sets for large multiplicities; implemented by `if any(spin_system.comp.mults>16,'all')`.
- Lines 64-65: Count chemical substances; implemented by `nsubst=numel(spin_system.chem.parts)`.
- Lines 67-113: Run connectivity analysis for IK-DNP basis set; builds `ee_conmatrix`, `en_conmatrix`, and `nn_conmatrix` from the coupling tensors above `spin_system.tols.inter_cutoff`.
- Lines 115-136: Coupling tensor norm for the connectivity analysis; `tensor_norm` is the isotropic part or the 2-norm of a coupling tensor according to `bas.connectivity`, shared by `IK-1`, `IK-2`, and `IK-SBS`.
- Lines 138-165: Run connectivity analysis for IK-1,2 basis sets; builds `spin_system.inter.conmatrix` from the spin-spin coupling tensors above `tols.inter_cutoff`, symmetrises it and `spin_system.inter.proxmatrix`, reports their densities and the number of uncoupled subsystems.
- Lines 167-278: Run connectivity analysis for IK-SBS basis set; requires spins and modes, builds the spin-spin graph and adds bosonic mode couplings above the cut-off (pairwise channels by 2-norm, `coupling_mod` and `zeeman_mod` inner blocks linked to their modes), symmetrises, and splits into `bb_conmatrix`, `sb_conmatrix`, and `ss_conmatrix` with diagonals set; reports the three densities.
- Lines 280-284: Build state lists for individual spins; implemented by `spin_state_lists{n}=(0:(spin_system.comp.mults(n)^2-1))'`.
- Lines 286-313: Apply longitudinal filters; for each substance, resolves `bas.longitudinal{s}` entries to spins of that substance and removes states with `M~=0` from their lists.
- Lines 315-316: Compute subspace dimensions for individual spins; implemented by `spin_dims=cellfun(@numel,spin_state_lists)`.
- Lines 318-319: Preallocate subgraph lists and their substance indices; implemented by `subgraphs=cell(nsubst,1); subgraph_subst=cell(nsubst,1);`
- Lines 321-481: Loop over chemical substances; generates coupling and proximity subgraphs in the spin index of the substance (`none`, `IK-0`, `IK-1`, `IK-2`, `IK-DNP`, and `IK-SBS` with `dfpt` on the three split graphs to `bas.inter_level(1:3)`), adds the rows of `bas.manual` that belong to the substance, removes spin-zero particles, empty, identical, and enclosed subgraphs, and embeds the result into the full spin index.
- Lines 483-485: Merge the subgraph lists of all substances; implemented by `subgraphs=vertcat(subgraphs{:})`.
- Lines 487-515: Resolve zero-quantum filter spins for each substance; implemented by `zq_spins(s,spins_in_question)=true()`.
- Lines 517-519: Balance the subgraph list; implemented by `randperm`.
- Lines 521-606: Populate the basis descriptor array; `parfor` over subgraphs builds the direct product descriptor, applies the coherence order, zero-quantum, IK-DNP inter-nuclear, and IK-SBS pure boson-boson and pure spin-spin correlation filters, and embeds it into the full spin index with `sparse`.
- Lines 608-609: Deallocate variables; implemented by `clear('spin_state_lists','subgraphs','subgraph_subst','spin_dims','zq_spins');`
- Lines 611-612: Pull basis descriptor from the nodes, unit state first; implemented by `basis_spec=[sparse(1,spin_system.comp.nspins); vertcat(basis_spec{:})];`
- Lines 614-616: Eliminate redundant states using a hash table; implemented by `basis_spec=unihash(basis_spec)`.
- Lines 618-632: Sort the basis explicitly; implemented by `sortrows`, distributed with `distrib_dim` for large bases.
- Lines 634-635: Deallocate variables; implemented by `clear('basis_spec');`
- Lines 637-641: Total projection quantum number and correlation order of each state; implemented by `[L,M]=lin2lm(spin_system.bas.basis)`.
- Lines 643-647: Report on chemical species; state counts per substance.
- Lines 649-650: Print the summary; implemented by `summary_basis(spin_system);`
- Lines 652-655: Run the symmetry treatment; implemented by `spin_system=symmetry(spin_system,bas);`
- Lines 657-682: Process Hilbert space Zeeman basis; Zeeman index table for `zeeman-hilb` and `zeeman-wavef`.
- Lines 660-661: Preallocate basis set array; implemented by `spin_system.bas.basis=zeros(prod(spin_system.comp.mults),spin_system.comp.nspins);`
- Lines 663-674: Fill basis set array; implemented by `for n=1:spin_system.comp.nspins`
- Lines 676-677: Report to the user; implemented by `report(spin_system,['matrix dimension for all operators and states: ' num2str(prod(spin_system.comp.`
- Lines 679-682: Run the symmetry treatment; implemented by `spin_system=symmetry(spin_system,bas);`
- Lines 684-711: Process Liouville space Zeeman basis; ket and bra index tables for `zeeman-liouv`.
- Lines 687-700: Build the Hilbert space Zeeman index table; implemented by `dim=prod(spin_system.comp.mults);`
- Lines 702-703: Ket and bra index tables in the vectorisation order; implemented by `spin_system.bas.basis=[repmat(zbas,[dim 1]) kron(zbas,ones(dim,1))];`
- Lines 705-706: Report to the user; implemented by `report(spin_system,['matrix dimension for all superoperators and state vectors: ' num2str(dim^2)]);`
- Lines 708-711: Run the symmetry treatment; implemented by `spin_system=symmetry(spin_system,bas);`
- Lines 713-740: Preload Lie algebra structure tables; implemented by `ist_product_table`.
- Lines 716-717: Inform the user; implemented by `report(spin_system,'caching Lie structure tables...');`
- Lines 719-720: Find the spin multiplicities present; implemented by `unique_mults=unique(spin_system.comp.mults);`
- Lines 722-724: Preallocate the structure table arrays; implemented by `spin_system.bas.lpst=cell(max(unique_mults),1);`
- Lines 726-740: Fill the arrays; implemented by `for n=setdiff(unique_mults,1)`
- Lines 742-748: Hash the basis descriptor for caching tools later; implemented by `md5_hash(spin_system.bas.basis)`.

### Control flow inferred from the code

- Line 57: conditional branch on `strcmp(spin_system.bas.formalism,'sphten-liouv')`.
- Line 68: conditional branch on `strcmp(spin_system.bas.approximation,'IK-DNP')`.
- Line 116: conditional branch on `ismember(spin_system.bas.approximation,{'IK-1','IK-2','IK-SBS'})` (coupling tensor norm).
- Line 139: conditional branch on `ismember(spin_system.bas.approximation,{'IK-1','IK-2'})` (spin-only coupling graph).
- Line 168: conditional branch on `strcmp(spin_system.bas.approximation,'IK-SBS')` (spin-boson coupling graphs).
- Line 287: conditional branch on `isfield(bas,'longitudinal')`.
- Line 322: `for` loop over `s=1:nsubst`.
- Line 330: dispatches on `spin_system.bas.approximation`; cases `'none'`, `'IK-0'`, `'IK-1'`, `'IK-2'`, `'IK-DNP'`, `'IK-SBS'`.
- Line 443: conditional branch on `isfield(bas,'manual')`; a row of `bas.manual` that spans two substances is an error.
- Line 489: conditional branch on `isfield(bas,'zero_quantum')`.
- Line 524: `parfor` loop over `n=1:size(subgraphs,1)`.
- Line 658: conditional branch on `ismember(spin_system.bas.formalism,{'zeeman-hilb','zeeman-wavef'})`.
- Line 685: conditional branch on `strcmp(spin_system.bas.formalism,'zeeman-liouv')`.

### Local helper functions

- Line 751: `grumble()` — `function grumble(spin_system,bas)`. Checks `bas.formalism`, `bas.approximation`, `bas.connectivity`, `bas.inter_level`, `bas.prox_level`, `bas.manual`, `bas.projections`, `bas.longitudinal`, and `bas.zero_quantum`; the three filters must be cell arrays with one element per chemical substance and are only available for `sphten-liouv`.

## Parameters / inputs

- spin_system - primary Spinach data structure, the output of create.m function
- bas - basis set specification structure described in detail in the online manual

## Outputs

- spin_system - primary Spinach data structure, updated with the basis set and related information
- Note: it is important to understand the factors that influence basis set selection in spin dynamics simulations - see our paper for further information on this subject.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `banner()`, `grumble()`, `cellfun()`, `summary_basis_opts()`, `report()`, `dfpt()`, `nchoosek()`, `prune_subgraphs()`, `unique()`, `lin2lm()`, `repelem()`, `sparse()`, `unihash()`, `distrib_dim()`, `sortrows()`, `summary_basis()`, `symmetry()`, `ist_product_table()`, `md5_hash()`.
