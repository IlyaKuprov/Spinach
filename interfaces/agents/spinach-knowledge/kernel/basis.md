# kernel/basis.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/basis.m`
- Signature: `spin_system=basis(spin_system,bas)`
- Total lines: 821

## Purpose

Basis set control. This is the second mandatory function (after create.m) that must be called in every calculation to build spin_system data struc- ture. Syntax: spin_system=basis(spin_system,bas)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 33-34: Show the banner; implemented by `banner(spin_system,'basis_banner')`.
- Lines 36-37: Check the input; implemented by `grumble(spin_system,bas)`.
- Lines 39-40: Store the settings; implemented by `spin_system.bas=bas`.
- Lines 42-43: Find electrons and nuclei; implemented by `e_idx=cellfun(@iselectron,spin_system.comp.isotopes)`.
- Lines 46-47: Report back to the user; implemented by `summary_basis_opts(spin_system)`.
- Lines 49-50: Process spherical tensor basis sets; implemented by `if strcmp(spin_system.bas.formalism,'sphten-liouv')`.
- Lines 52-53: Disallow spherical tensor basis sets for large multiplicities; implemented by `if any(spin_system.comp.mults>16,'all')`.
- Lines 57-58: Run connectivity analysis for IK-DNP basis set; implemented by `if strcmp(spin_system.bas.approximation,'IK-DNP')`.
- Lines 60-61: Make sure there are only electrons and nuclei; implemented by `if (nnz(e_idx)==0)||(nnz(n_idx)==0)`.
- Lines 68-69: Isolate three types of interactions (e-e, n-n, e-n); implemented by `ee_couplings=spin_system.inter.coupling.matrix`.
- Lines 76-78: Remind the user about the amplitude cut-off; implemented by `report(spin_system,['coupling tensors with norm below ' num2str(spin_system.tols.inter_cutoff) ' Hz will be ignored.'])`.
- Lines 80-81: Generate three types of connectivity graphs (e-e, e-n, n-n); implemented by `ee_conmatrix=sparse(cellfun(@(x)norm(x,2),ee_couplings)>2*pi*spin_system.tols.inter_cutoff)`.
- Lines 85-86: Make sure each spin is connected to itself; implemented by `ee_conmatrix=ee_conmatrix|speye(size(ee_conmatrix))`.
- Lines 90-91: Make sure connectivity is reciprocal; implemented by `ee_conmatrix=ee_conmatrix|transpose(ee_conmatrix)`.
- Lines 97-98: Run connectivity analysis for IK-1,2 basis sets; implemented by `if ismember(spin_system.bas.approximation,{'IK-1','IK-2'})`.
- Lines 100-101: Build connectivity and proximity matrices; implemented by `switch bas.connectivity`.
- Lines 105-106: Use scalar parts of all interaction tensors; implemented by `report(spin_system,'scalar couplings will be used to build the coupling graph.')`.
- Lines 111-112: Use complete interaction tensors; implemented by `report(spin_system,'full coupling tensors will be used to build the coupling graph.')`.

### Control flow inferred from the code

- Line 50: conditional branch on `strcmp(spin_system.bas.formalism,'sphten-liouv')`.
- Line 53: conditional branch on `any(spin_system.comp.mults>16,'all')`.
- Line 58: conditional branch on `strcmp(spin_system.bas.approximation,'IK-DNP')`.
- Line 61: conditional branch on `(nnz(e_idx)==0)||(nnz(n_idx)==0)`.
- Line 64: conditional branch on `~all(e_idx|n_idx)`.
- Line 98: conditional branch on `ismember(spin_system.bas.approximation,{'IK-1','IK-2'})`.
- Line 101: dispatches on `bas.connectivity`; cases `'scalar_couplings'`, `'full_tensors'`, `'none'`, `'IK-0'`, `'IK-1'`, `'IK-2'`, `'IK-DNP'`.
- Line 118: conditional branch on `isfield(spin_system.inter,'modes')`.
- Line 123: `for` loop over `n=1:numel(chan_flds)`.
- Line 148: conditional branch on `n_subsystems>1`.
- Line 156: `for` loop over `n=1:spin_system.comp.nspins`.
- Line 161: conditional branch on `isfield(bas,'longitudinals')`.
- Line 165: `for` loop over `n=1:numel(bas.longitudinals)`.
- Line 166: conditional branch on `isnumeric(bas.longitudinals{n})`.

### Key state/data transformations

- Lines 40: computes `spin_system.bas` using `spin_system.bas=bas`.
- Lines 43: computes `e_idx` using `e_idx=cellfun(@iselectron,spin_system.comp.isotopes)`.
- Lines 44: computes `n_idx` using `n_idx=cellfun(@isnucleus,spin_system.comp.isotopes)`.
- Lines 69: computes `ee_couplings` using `ee_couplings=spin_system.inter.coupling.matrix`.
- Lines 70: computes `ee_couplings(:,n_idx)` using `ee_couplings(:,n_idx)={[]}; ee_couplings(n_idx,:)={[]}`.
- Lines 71: computes `en_couplings` using `en_couplings=spin_system.inter.coupling.matrix`.
- Lines 72: computes `en_couplings(e_idx,e_idx)` using `en_couplings(e_idx,e_idx)={[]}; en_couplings(n_idx,n_idx)={[]}`.
- Lines 73: computes `nn_couplings` using `nn_couplings=spin_system.inter.coupling.matrix`.
- Lines 74: computes `nn_couplings(:,e_idx)` using `nn_couplings(:,e_idx)={[]}; nn_couplings(e_idx,:)={[]}`.
- Lines 81: computes `ee_conmatrix` using `ee_conmatrix=sparse(cellfun(@(x)norm(x,2),ee_couplings)>2*pi*spin_system.tols.inter_cutoff)`.
- Lines 82: computes `en_conmatrix` using `en_conmatrix=sparse(cellfun(@(x)norm(x,2),en_couplings)>2*pi*spin_system.tols.inter_cutoff)`.
- Lines 83: computes `nn_conmatrix` using `nn_conmatrix=sparse(cellfun(@(x)norm(x,2),nn_couplings)>2*pi*spin_system.tols.inter_cutoff)`.
- Lines 107: computes `spin_system.inter.conmatrix` using `spin_system.inter.conmatrix=sparse(abs(cellfun(@trace,spin_system.inter.coupling.matrix)/3)>2*pi*spin_system.tols.inter_cutoff)`.
- Lines 120-121: computes `chan_flds` using `chan_flds={'exchange','dispersive','kerr','longitudinal', 'coupling_mod','zeeman_mod'}`.
- Lines 122: computes `mode_conmat` using `mode_conmat=false(spin_system.comp.nspins)`.
- Lines 134: computes `spin_system.inter.proxmatrix` using `spin_system.inter.proxmatrix=spin_system.inter.proxmatrix|speye(size(spin_system.inter.proxmatrix))`.
- Lines 145: computes `n_subsystems` using `n_subsystems=max(scomponents(spin_system.inter.conmatrix|spin_system.inter.proxmatrix))`.
- Lines 155: computes `spin_state_lists` using `spin_state_lists=cell(spin_system.comp.nspins,1)`.

### Local helper functions

- Line 610: `grumble()` — `function grumble(spin_system,bas)`. Check bas.formalism
  - Representative operation: `if ~isfield(bas,'formalism')`.
  - Representative operation: `error('basis specification in bas.formalism is required.')`.

## Parameters / inputs

- spin_system -primary Spinach data structure, the output
- of create.m function
- bas -basis set specification structure described
- in detail in the online manual

## Outputs

- spin_system -primary Spinach data structure, updated with
- the basis set and related information
- Note: it is important to understand the factors that influence basis set
- selection in spin dynamics simulations -see our paper
- for further information on this subject.

## Implementation structure

- Basis set control. This is the second mandatory function (after create.m)
- that must be called in every calculation to build spin_system data struc-
- ture. Syntax:
- spin_system=basis(spin_system,bas)
- spin_system -primary Spinach data structure, the output
- of create.m function
- bas -basis set specification structure described
- in detail in the online manual
- spin_system -primary Spinach data structure, updated with
- the basis set and related information
- Note: it is important to understand the factors that influence basis set
- selection in spin dynamics simulations -see our paper

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `banner()`, `grumble()`, `cellfun()`, `summary_basis_opts()`, `strcmp()`, `any()`, `nnz()`, `all()`, `ee_couplings()`, `en_couplings()`, `nn_couplings()`, `report()`, `num2str()`, `speye()`, `transpose()`, `ismember()`.
