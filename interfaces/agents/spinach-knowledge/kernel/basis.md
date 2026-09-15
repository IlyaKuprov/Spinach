# kernel/basis.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/basis.m`
- Signature: `spin_system=basis(spin_system,bas)`
- Total lines: 906

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

- Lines 32-33: Show the banner; implemented by `banner(spin_system,'basis_banner')`.
- Lines 35-36: Check the input; implemented by `grumble(spin_system,bas)`.
- Lines 38-39: Store the settings; implemented by `spin_system.bas=bas`.
- Lines 41-43: Find electrons and nuclei; implemented by `e_idx=cellfun(@iselectron,spin_system.comp.isotopes)`.
- Lines 45-46: Report back to the user; implemented by `summary_basis_opts(spin_system)`.
- Lines 48-51: Remind the user about the amplitude cut-off; implemented by `report(spin_system,['coupling tensors with norm below ' ...`
- Lines 53-565: Process spherical tensor basis sets; implemented by `if strcmp(spin_system.bas.formalism,'sphten-liouv')`.
- Lines 56-59: Disallow spherical tensor basis sets for large multiplicities; implemented by `if any(spin_system.comp.mults>16,'all')`.
- Lines 61-62: Count chemical substances; implemented by `nsubst=numel(spin_system.chem.parts)`.
- Lines 64-110: Run connectivity analysis for IK-DNP basis set; builds `ee_conmatrix`, `en_conmatrix`, and `nn_conmatrix` from the coupling tensors above `spin_system.tols.inter_cutoff`.
- Lines 112-231: Run connectivity analysis for IK-1,2 basis sets; builds `spin_system.inter.conmatrix` from the coupling tensors above `tols.inter_cutoff` under the `bas.connectivity` norm, adds bosonic mode couplings above the same cut-off (pairwise `exchange`, `dispersive`, `kerr`, and `longitudinal` channels by 2-norm; `coupling_mod` and `zeeman_mod` entries connect the modes of the outer cell to the spin pairs and spins in the inner block whose derivative is above the cut-off), symmetrises it and `spin_system.inter.proxmatrix`, and reports their densities.
- Lines 233-237: Build state lists for individual spins; implemented by `spin_state_lists{n}=(0:(spin_system.comp.mults(n)^2-1))'`.
- Lines 239-266: Apply longitudinal filters; for each substance, resolves `bas.longitudinal{s}` entries to spins of that substance and removes states with `M~=0` from their lists.
- Lines 268-269: Compute subspace dimensions for individual spins; implemented by `spin_dims=cellfun(@numel,spin_state_lists)`.
- Lines 271-272: Preallocate subgraph lists and their substance indices; implemented by `subgraphs=cell(nsubst,1); subgraph_subst=cell(nsubst,1);`
- Lines 274-411: Loop over chemical substances; generates coupling and proximity subgraphs in the spin index of the substance (`none`, `IK-0`, `IK-1`, `IK-2`, `IK-DNP`), adds the rows of `bas.manual` that belong to the substance, removes spin-zero particles, empty, identical, and enclosed subgraphs, and embeds the result into the full spin index.
- Lines 413-415: Merge the subgraph lists of all substances; implemented by `subgraphs=vertcat(subgraphs{:})`.
- Lines 417-445: Resolve zero-quantum filter spins for each substance; implemented by `zq_spins(s,spins_in_question)=true()`.
- Lines 447-449: Balance the subgraph list; implemented by `randperm`.
- Lines 451-516: Populate the basis descriptor array; `parfor` over subgraphs builds the direct product descriptor, applies the coherence order, zero-quantum, and IK-DNP inter-nuclear filters, and embeds it into the full spin index with `sparse`.
- Lines 518-519: Deallocate variables; implemented by `clear('spin_state_lists','subgraphs','subgraph_subst','spin_dims','zq_spins');`
- Lines 521-522: Pull basis descriptor from the nodes, unit state first; implemented by `basis_spec=[sparse(1,spin_system.comp.nspins); vertcat(basis_spec{:})];`
- Lines 524-526: Eliminate redundant states using a hash table; implemented by `basis_spec=unihash(basis_spec)`.
- Lines 528-542: Sort the basis explicitly; implemented by `sortrows`, distributed with `distrib_dim` for large bases.
- Lines 544-545: Deallocate variables; implemented by `clear('basis_spec');`
- Lines 547-551: Total projection quantum number and correlation order of each state; implemented by `[L,M]=lin2lm(spin_system.bas.basis)`.
- Lines 553-557: Report on chemical species; state counts per substance.
- Lines 559-560: Print the summary; implemented by `summary_basis(spin_system);`
- Lines 562-565: Run the symmetry treatment; implemented by `spin_system=symmetry(spin_system,bas);`
- Lines 567-592: Process Hilbert space Zeeman basis; Zeeman index table for `zeeman-hilb` and `zeeman-wavef`.
- Lines 570-571: Preallocate basis set array; implemented by `spin_system.bas.basis=zeros(prod(spin_system.comp.mults),spin_system.comp.nspins);`
- Lines 573-584: Fill basis set array; implemented by `for n=1:spin_system.comp.nspins`
- Lines 586-587: Report to the user; implemented by `report(spin_system,['matrix dimension for all operators and states: ' num2str(prod(spin_sy`
- Lines 589-592: Run the symmetry treatment; implemented by `spin_system=symmetry(spin_system,bas);`
- Lines 594-621: Process Liouville space Zeeman basis; ket and bra index tables for `zeeman-liouv`.
- Lines 597-610: Build the Hilbert space Zeeman index table; implemented by `dim=prod(spin_system.comp.mults);`
- Lines 612-613: Ket and bra index tables in the vectorisation order; implemented by `spin_system.bas.basis=[repmat(zbas,[dim 1]) kron(zbas,ones(dim,1))];`
- Lines 615-616: Report to the user; implemented by `report(spin_system,['matrix dimension for all superoperators and state vectors: ' num2str(`
- Lines 618-621: Run the symmetry treatment; implemented by `spin_system=symmetry(spin_system,bas);`
- Lines 623-650: Preload Lie algebra structure tables; implemented by `ist_product_table`.
- Lines 626-627: Inform the user; implemented by `report(spin_system,'caching Lie structure tables...');`
- Lines 629-630: Find the spin multiplicities present; implemented by `unique_mults=unique(spin_system.comp.mults);`
- Lines 632-634: Preallocate the structure table arrays; implemented by `spin_system.bas.lpst=cell(max(unique_mults),1);`
- Lines 636-650: Fill the arrays; implemented by `for n=setdiff(unique_mults,1)`
- Lines 652-658: Hash the basis descriptor for caching tools later; implemented by `md5_hash(spin_system.bas.basis)`.

### Control flow inferred from the code

- Line 54: conditional branch on `strcmp(spin_system.bas.formalism,'sphten-liouv')`.
- Line 65: conditional branch on `strcmp(spin_system.bas.approximation,'IK-DNP')`.
- Line 113: conditional branch on `ismember(spin_system.bas.approximation,{'IK-1','IK-2'})`.
- Line 240: conditional branch on `isfield(bas,'longitudinal')`.
- Line 275: `for` loop over `s=1:nsubst`.
- Line 283: dispatches on `spin_system.bas.approximation`; cases `'none'`, `'IK-0'`, `'IK-1'`, `'IK-2'`, `'IK-DNP'`.
- Line 373: conditional branch on `isfield(bas,'manual')`; a row of `bas.manual` that spans two substances is an error.
- Line 419: conditional branch on `isfield(bas,'zero_quantum')`.
- Line 454: `parfor` loop over `n=1:size(subgraphs,1)`.
- Line 568: conditional branch on `ismember(spin_system.bas.formalism,{'zeeman-hilb','zeeman-wavef'})`.
- Line 595: conditional branch on `strcmp(spin_system.bas.formalism,'zeeman-liouv')`.

### Local helper functions

- Line 661: `grumble()` — `function grumble(spin_system,bas)`. Checks `bas.formalism`, `bas.approximation`, `bas.connectivity`, `bas.inter_level`, `bas.prox_level`, `bas.manual`, `bas.projections`, `bas.longitudinal`, and `bas.zero_quantum`; the three filters must be cell arrays with one element per chemical substance and are only available for `sphten-liouv`.

## Parameters / inputs

- spin_system - primary Spinach data structure, the output of create.m function
- bas - basis set specification structure described in detail in the online manual

## Outputs

- spin_system - primary Spinach data structure, updated with the basis set and related information
- Note: it is important to understand the factors that influence basis set selection in spin dynamics simulations - see our paper for further information on this subject.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `banner()`, `grumble()`, `cellfun()`, `summary_basis_opts()`, `report()`, `dfpt()`, `nchoosek()`, `prune_subgraphs()`, `unique()`, `lin2lm()`, `repelem()`, `sparse()`, `unihash()`, `distrib_dim()`, `sortrows()`, `summary_basis()`, `symmetry()`, `ist_product_table()`, `md5_hash()`.
