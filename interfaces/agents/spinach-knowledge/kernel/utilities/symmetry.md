# kernel/utilities/symmetry.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/symmetry.m`
- Signature: `spin_system=symmetry(spin_system,bas)`
- Total lines: 351

## Purpose

Permutation symmetry treatment. Compiles character tables of composite symmetry groups, builds the permutation table for each spin state in the basis, and builds projectors into the irreducible representations of the direct product symmetry group. Syntax: spin_system=symmetry(spin_system,bas)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 42-43: Check consistency; implemented by `grumble(spin_system,bas)`.
- Lines 45-46: Check the disable switch; implemented by `if ismember('symmetry',spin_system.sys.disable)`.
- Lines 48-49: Issue a reminder to the user; implemented by `report(spin_system,'WARNING - symmetry factorization disabled by the user.')`.
- Lines 51-52: Write empty cells; implemented by `spin_system.comp.sym_group={}`.
- Lines 58-59: Symmetry group; implemented by `if isfield(bas,'sym_group')`.
- Lines 65-66: Symmetry-related spins; implemented by `if isfield(bas,'sym_spins')`.
- Lines 72-73: Irreducible representation composition; implemented by `if isfield(bas,'sym_a1g_only')`.
- Lines 81-82: Report back to the user; implemented by `if ~isempty(spin_system.comp.sym_group)`.
- Lines 86-87: Validate that the interactions respect the declared symmetry; implemented by `validate_sym(spin_system,bas)`.
- Lines 89-90: Compute group direct product if necessary; implemented by `if numel(spin_system.comp.sym_group)>1`.
- Lines 92-93: Lift constituent groups from the database; implemented by `ngroups=numel(spin_system.comp.sym_group); groups=cell(1,ngroups)`.
- Lines 98-99: Compute direct product character table; implemented by `group.characters=1`.
- Lines 108-109: Compute direct product element list; implemented by `group.elements=groups{1}.elements; group.order=groups{1}.order`.
- Lines 117-118: Concatenate spin lists; implemented by `spins=horzcat(spin_system.comp.sym_spins{:})`.
- Lines 122-123: Lift the group from the database; implemented by `spins=spin_system.comp.sym_spins{1}`.
- Lines 130-131: Remind the user that symmetry is not operational; implemented by `report(spin_system,'no symmetry information available.')`.
- Lines 135-136: Run the SALC procedure; implemented by `if exist('group','var')`.
- Lines 138-139: Preallocate the permutation table; implemented by `permutation_table=zeros(size(spin_system.bas.basis,1),group.order)`.

### Control flow inferred from the code

- Line 46: conditional branch on `ismember('symmetry',spin_system.sys.disable)`.
- Line 59: conditional branch on `isfield(bas,'sym_group')`.
- Line 66: conditional branch on `isfield(bas,'sym_spins')`.
- Line 73: conditional branch on `isfield(bas,'sym_a1g_only')`.
- Line 82: conditional branch on `~isempty(spin_system.comp.sym_group)`.
- Line 90: conditional branch on `numel(spin_system.comp.sym_group)>1`.
- Line 94: `for` loop over `n=1:ngroups`.
- Line 100: `for` loop over `n=1:ngroups`.
- Line 110: `for` loop over `n=2:ngroups`.
- Line 136: conditional branch on `exist('group','var')`.
- Line 142: `parfor` loop over `n=1:group.order`.
- Line 145: conditional branch on `strcmp(spin_system.bas.formalism,'zeeman-liouv')`.
- Line 154: conditional branch on `spin_system.comp.sym_a1g_only`.
- Line 188: `for` loop over `n=1:group.n_irreps`.

### Key state/data transformations

- Lines 52: computes `spin_system.comp.sym_group` using `spin_system.comp.sym_group={}`.
- Lines 53: computes `spin_system.comp.sym_spins` using `spin_system.comp.sym_spins={}`.
- Lines 54: computes `spin_system.comp.sym_a1g_only` using `spin_system.comp.sym_a1g_only=true()`.
- Lines 93: computes `ngroups` using `ngroups=numel(spin_system.comp.sym_group); groups=cell(1,ngroups)`.
- Lines 95: computes `groups{n}` using `groups{n}=perm_group(spin_system.comp.sym_group{n})`.
- Lines 99: computes `group.characters` using `group.characters=1`.
- Lines 103: computes `group.irrep_dims` using `group.irrep_dims=group.characters(:,1)'`.
- Lines 104: computes `group.n_irreps` using `group.n_irreps=size(group.characters,1)`.
- Lines 109: computes `group.elements` using `group.elements=groups{1}.elements; group.order=groups{1}.order`.
- Lines 113: computes `group.order` using `group.order=group.order*groups{n}.order`.
- Lines 118: computes `spins` using `spins=horzcat(spin_system.comp.sym_spins{:})`.
- Lines 124: computes `group` using `group=perm_group(spin_system.comp.sym_group{1})`.
- Lines 139: computes `permutation_table` using `permutation_table=zeros(size(spin_system.bas.basis,1),group.order)`.
- Lines 143: computes `group_element` using `group_element=1:spin_system.comp.nspins`.
- Lines 144: computes `group_element(spins)` using `group_element(spins)=group_element(spins(group.elements(n,:)))`.
- Lines 148: computes `permuted_basis` using `permuted_basis=spin_system.bas.basis(:,group_element)`.
- Lines 149: computes `index` using `index=spsortrows(sparse(permuted_basis))`.
- Lines 150: computes `permutation_table(:,n)` using `permutation_table(:,n)=index`.

### Local helper functions

- Line 276: `grumble()` — `function grumble(spin_system,bas)`. Check symmetry parameters
  - Representative operation: `if isfield(bas,'sym_group')`.
  - Representative operation: `if (~iscell(bas.sym_group))||any(~cellfun(@ischar,bas.sym_group))`.

## Parameters / inputs

- spin_system -Spinach spin system description object
- produced as described in the spin system
- and basis specification sections, of the
- of the online manual.
- bas -basis input structure described in the
- basis specification section of the manual

## Outputs

- spin_system.bas.irrep(n).projector -projector matrices
- into each irreducible
- representation
- spin_system.bas.irrep(n).dimension -dimension of each ir-
- reducible representa-
- tion
- Note: this is a service function of the Spinach kernel that
- should not be called directly; it is called by basis.m
- Note: non-Abelian groups and multi-dimensional irreps are sup-
- ported -edit perm_group.m to add your own groups.

## Implementation structure

- Permutation symmetry treatment. Compiles character tables of
- composite symmetry groups, builds the permutation table for
- each spin state in the basis, and builds projectors into the
- irreducible representations of the direct product symmetry
- group. Syntax:
- spin_system=symmetry(spin_system,bas)
- spin_system - Spinach spin system description object
- produced as described in the spin system
- and basis specification sections, of the
- of the online manual.
- bas - basis input structure described in the
- basis specification section of the manual

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ismember()`, `report()`, `true()`, `isfield()`, `false()`, `summary_symmetry()`, `validate_sym()`, `perm_group()`, `num2str()`, `horzcat()`, `isscalar()`, `exist()`, `group_element()`, `spins()`, `strcmp()`.
