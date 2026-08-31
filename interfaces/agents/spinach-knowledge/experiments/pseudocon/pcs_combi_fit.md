# experiments/pseudocon/pcs_combi_fit.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/pseudocon/pcs_combi_fit.m`
- Signature: `[d_shifts,p_shifts,pcs_theo,pcs_expt,chi,total_theo]=pcs_combi_fit(parameters)`
- Total lines: 187

## Purpose

Combinatorial PCS fitting function. Takes into account potential am- biguities in diamagnetic and paramagnetic NMR assignments. Syntax: [d_shifts,p_shifts,pcs_theo,... pcs_expt,chi,total_theo]=pcs_combi_fit(parameters)

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 62-63: Check consistency; implemented by `grumble(parameters)`.
- Lines 65-66: Build diamagnetic shift permutation table; implemented by `perm_table_d=perms(parameters.d_ambig{1})`.
- Lines 78-79: Build paramagnetic shift permutation table; implemented by `perm_table_p=perms(parameters.p_ambig{1})`.
- Lines 91-92: Compute permutation scores; implemented by `N=size(diamag_shift_perms,1); M=size(paramag_shift_perms,1); score=zeros(N*M,1)`.
- Lines 96-97: Disentangle indices; implemented by `[d_index,p_index]=ind2sub([N M],n)`.
- Lines 99-100: Compute pseudocontact shifts; implemented by `pcs=paramag_shift_perms(p_index,:)-diamag_shift_perms(d_index,:)`.
- Lines 102-103: Build the input arrays; implemented by `hfcs=cell(sum(cellfun(@numel,parameters.spin_groups(:))),1); p=1`.
- Lines 112-113: Get the assignment score; implemented by `[~,score(n)]=pcs2chi(hfcs,pcs_expt,parameters.isotopes)`.
- Lines 117-118: Find the shift permutations that led to the best score; implemented by `[min_score,n]=min(score)`.
- Lines 134-135: Recover the susceptibility; implemented by `chi=pcs2chi(hfcs,pcs_expt,parameters.isotopes)`.
- Lines 137-138: Compute the predicted pseudocontact shifts; implemented by `pcs_theo=zeros(numel(hfcs),1)`.
- Lines 143-144: Compute the predicted total shifts; implemented by `total_theo=pcs_theo+repelem(d_shifts(:),cellfun(@numel,parameters.spin_groups(:)))`.
- Lines 146-147: Report to the user; implemented by `disp('Best permuted diamagnetic shifts:'); disp(d_shifts')`.

### Control flow inferred from the code

- Line 67: `for` loop over `n=2:numel(parameters.d_ambig)`.
- Line 73: `for` loop over `n=1:size(perm_table_d,1)`.
- Line 80: `for` loop over `n=2:numel(parameters.p_ambig)`.
- Line 86: `for` loop over `n=1:size(perm_table_p,1)`.
- Line 94: `parfor` loop over `n=1:N*M`.
- Line 105: `for` loop over `k=1:numel(parameters.spin_groups)`.
- Line 106: `for` loop over `m=1:numel(parameters.spin_groups{k})`.
- Line 127: `for` loop over `k=1:numel(parameters.spin_groups)`.
- Line 128: `for` loop over `m=1:numel(parameters.spin_groups{k})`.
- Line 139: `for` loop over `k=1:numel(hfcs)`.

### Key state/data transformations

- Lines 66: computes `perm_table_d` using `perm_table_d=perms(parameters.d_ambig{1})`.
- Lines 68: computes `local_perms` using `local_perms=perms(parameters.d_ambig{n})`.
- Lines 72: computes `diamag_shift_perms` using `diamag_shift_perms=zeros(size(perm_table_d,1),numel(parameters.d_shifts))`.
- Lines 74: computes `diamag_shift_perms(n,:)` using `diamag_shift_perms(n,:)=parameters.d_shifts`.
- Lines 75: computes `diamag_shift_perms(n,cell2mat(parameters.d_ambig))` using `diamag_shift_perms(n,cell2mat(parameters.d_ambig))=diamag_shift_perms(n,perm_table_d(n,:))`.
- Lines 79: computes `perm_table_p` using `perm_table_p=perms(parameters.p_ambig{1})`.
- Lines 85: computes `paramag_shift_perms` using `paramag_shift_perms=zeros(size(perm_table_p,1),numel(parameters.p_shifts))`.
- Lines 87: computes `paramag_shift_perms(n,:)` using `paramag_shift_perms(n,:)=parameters.p_shifts`.
- Lines 88: computes `paramag_shift_perms(n,cell2mat(parameters.p_ambig))` using `paramag_shift_perms(n,cell2mat(parameters.p_ambig))=paramag_shift_perms(n,perm_table_p(n,:))`.
- Lines 92: computes `N` using `N=size(diamag_shift_perms,1); M=size(paramag_shift_perms,1); score=zeros(N*M,1)`.
- Lines 97: computes `[d_index,p_index]` using `[d_index,p_index]=ind2sub([N M],n)`.
- Lines 100: computes `pcs` using `pcs=paramag_shift_perms(p_index,:)-diamag_shift_perms(d_index,:)`.
- Lines 103: computes `hfcs` using `hfcs=cell(sum(cellfun(@numel,parameters.spin_groups(:))),1); p=1`.
- Lines 104: computes `pcs_expt` using `pcs_expt=zeros(sum(cellfun(@numel,parameters.spin_groups(:))),1)`.
- Lines 107: computes `hfcs{p}` using `hfcs{p}=parameters.hfcs{parameters.spin_groups{k}(m)}`.
- Lines 108: computes `pcs_expt(p)` using `pcs_expt(p)=pcs(k); p=p+1`.
- Lines 113: computes `[~,score(n)]` using `[~,score(n)]=pcs2chi(hfcs,pcs_expt,parameters.isotopes)`.
- Lines 118: computes `[min_score,n]` using `[min_score,n]=min(score)`.

### Local helper functions

- Line 156: `grumble()` — `function grumble(parameters)`.
  - Representative operation: `if ~isfield(parameters,'hfcs')`.
  - Representative operation: `error('parameters.hfcs field is missing.')`.

## Parameters / inputs

- parameters.hfcs -cell array of 3x3 hyperfine tensors,
- in Gauss, usually out of gparse() or
- something similar
- parameters.isotopes -cell array of isotope specificati-
- ons, e.g. {'1H','1H'}
- parameters.spin_groups -cell array of integer vectors
- specifying the numbers of spins
- that have each of the chemical
- shifts specified, e.g.
- {[28 22 30]; [25 33 81]}
- parameters.d_shifts -a vector of unique diamagnetic che-
- mical shifts, in ppm
- parameters.p_shifts -a vector of unique paramagnetic
- chemical shifts, in ppm
- parameters.d_ambig -a cell array of integer vectors spe-
- cifying the spins for which the dia-
- magnetic assignment can potentially
- be swapped around.
- parameters.p_ambig -a cell array of integer vectors spe-
- cifying the spins for which the para-
- magnetic assignment can potentially
- be swapped around.

## Outputs

- d_shifts -diamagnetic chemical shifts, optimally permuted
- p_shifts -paramagnetic chemical shifts, optimally permuted
- pcs_theo -theoretical pseudocontact shifts
- pcs_expt -experimental pseudocontact shifts from optimally
- permuted assignments
- chi -rank 2 part of the magnetic susceptibility tensor,
- in cubic Angstrom
- total_theo -theoretical total NMR chemcial shifts, computed
- as a sum of d_shifts and pcs_theo

## Implementation structure

- Combinatorial PCS fitting function. Takes into account potential am-
- biguities in diamagnetic and paramagnetic NMR assignments. Syntax:
- [d_shifts,p_shifts,pcs_theo,...
- pcs_expt,chi,total_theo]=pcs_combi_fit(parameters)
- parameters.hfcs -cell array of 3x3 hyperfine tensors,
- in Gauss, usually out of gparse() or
- something similar
- parameters.isotopes -cell array of isotope specificati-
- ons, e.g. {'1H','1H'}
- parameters.spin_groups -cell array of integer vectors
- specifying the numbers of spins
- that have each of the chemical

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `perms()`, `num2str()`, `diamag_shift_perms()`, `cell2mat()`, `perm_table_d()`, `paramag_shift_perms()`, `perm_table_p()`, `ind2sub()`, `cellfun()`, `pcs_expt()`, `pcs()`, `score()`, `pcs2chi()`, `p_shifts()`, `d_shifts()`.
