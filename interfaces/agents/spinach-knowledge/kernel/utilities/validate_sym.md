# kernel/utilities/validate_sym.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/validate_sym.m`
- Signature: `validate_sym(spin_system,bas)`
- Total lines: 149

## Purpose

Extended validation of user-declared permutation symmetry. Confirms that the Zeeman, coupling, and giant-spin interactions stored in the spin system object are strictly invariant under every operation of each declared permutation group, so that the irreducible representa- tion projectors built by symmetry.m correspond to a symmetry that the interactions actually possess. The declared symmetry is a permu- tation of sp

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 39-40: Check consistency; implemented by `grumble(spin_system,bas)`.
- Lines 42-43: Skip validation when symmetry is switched off; implemented by `if ismember('symmetry',spin_system.sys.disable), return; end`.
- Lines 45-46: Skip validation when no symmetry group is declared; implemented by `if (~isfield(bas,'sym_group'))||isempty(bas.sym_group), return; end`.
- Lines 48-49: Interaction agreement tolerance, rad/s; implemented by `tol=2*pi*spin_system.tols.inter_cutoff`.
- Lines 51-52: Number of spins in the system; implemented by `nspins=spin_system.comp.nspins`.
- Lines 54-55: Zeeman and giant-spin interaction arrays; implemented by `zeeman=spin_system.inter.zeeman.matrix`.
- Lines 58-59: Loop over declared symmetry groups; implemented by `for m=1:numel(bas.sym_group)`.
- Lines 61-62: Spins in the current symmetry group; implemented by `spins=bas.sym_spins{m}`.
- Lines 64-65: All spins in a symmetry group must be the same isotope; implemented by `if any(~strcmp(spin_system.comp.isotopes(spins),spin_system.comp.isotopes{spins(1)}))`.
- Lines 70-71: Permutation elements of the declared group; implemented by `group=perm_group(bas.sym_group{m})`.
- Lines 73-74: Loop over symmetry operations; implemented by `for n=1:group.order`.
- Lines 76-77: Global spin permutation for this operation; implemented by `perm=1:nspins; perm(spins)=spins(group.elements(n,:))`.
- Lines 79-80: Loop over spins in the symmetry group; implemented by `for k=1:numel(spins)`.
- Lines 82-83: Check Zeeman tensor invariance; implemented by `if norm(zeeman{perm(spins(k))}-zeeman{spins(k)},2)>tol`.
- Lines 89-90: Check giant-spin coefficient invariance across all ranks; implemented by `coeff_a=giant{spins(k)}; coeff_b=giant{perm(spins(k))}`.
- Lines 103-104: Check coupling tensor invariance against every spin; implemented by `for q=1:nspins`.
- Lines 120-121: Confirm success to the user; implemented by `report(spin_system,'declared permutation symmetry is consistent with the interaction data.')`.

### Control flow inferred from the code

- Line 43: conditional branch on `ismember('symmetry',spin_system.sys.disable), return; end`.
- Line 46: conditional branch on `(~isfield(bas,'sym_group'))||isempty(bas.sym_group), return; end`.
- Line 59: `for` loop over `m=1:numel(bas.sym_group)`.
- Line 65: conditional branch on `any(~strcmp(spin_system.comp.isotopes(spins),spin_system.comp.isotopes{spins(1)}))`.
- Line 74: `for` loop over `n=1:group.order`.
- Line 80: `for` loop over `k=1:numel(spins)`.
- Line 83: conditional branch on `norm(zeeman{perm(spins(k))}-zeeman{spins(k)},2)>tol`.
- Line 92: `for` loop over `r=1:numel(coeff_a)`.
- Line 93: conditional branch on `(~giant_mismatch)&&(norm(coeff_b{r}-coeff_a{r},2)>tol)`.
- Line 97: conditional branch on `giant_mismatch`.
- Line 104: `for` loop over `q=1:nspins`.
- Line 105: conditional branch on `norm(get_coupling(spin_system,perm(spins(k)),perm(q))-`.

### Key state/data transformations

- Lines 49: computes `tol` using `tol=2*pi*spin_system.tols.inter_cutoff`.
- Lines 52: computes `nspins` using `nspins=spin_system.comp.nspins`.
- Lines 55: computes `zeeman` using `zeeman=spin_system.inter.zeeman.matrix`.
- Lines 56: computes `giant` using `giant=spin_system.inter.giant.coeff`.
- Lines 62: computes `spins` using `spins=bas.sym_spins{m}`.
- Lines 71: computes `group` using `group=perm_group(bas.sym_group{m})`.
- Lines 77: computes `perm` using `perm=1:nspins; perm(spins)=spins(group.elements(n,:))`.
- Lines 90: computes `coeff_a` using `coeff_a=giant{spins(k)}; coeff_b=giant{perm(spins(k))}`.
- Lines 94: computes `giant_mismatch` using `giant_mismatch=true`.

### Local helper functions

- Line 126: `grumble()` — `function grumble(spin_system,bas)`.
  - Representative operation: `if isfield(bas,'sym_group')`.
  - Representative operation: `if (~iscell(bas.sym_group))||any(~cellfun(@ischar,bas.sym_group))`.

## Parameters / inputs

- spin_system -Spinach spin system description object, with the
- interaction arrays already processed by create.m
- bas -basis specification structure; the fields used
- here are bas.sym_group (cell array of group names)
- and bas.sym_spins (cell array of spin index vec-
- tors), as described in the basis specification
- section of the online manual

## Outputs

- this function returns nothing; it throws a descriptive error when
- the interaction data does not obey a declared symmetry
- Note: this is a service function of the Spinach kernel that should
- not be called directly; it is called by symmetry.m

## Implementation structure

- Extended validation of user-declared permutation symmetry. Confirms
- that the Zeeman, coupling, and giant-spin interactions stored in the
- spin system object are strictly invariant under every operation of
- each declared permutation group, so that the irreducible representa-
- tion projectors built by symmetry.m correspond to a symmetry that
- the interactions actually possess. The declared symmetry is a permu-
- tation of spin labels, not a spatial rotation; interaction tensors
- related by a rotation rather than being identical are not accepted.
- When no symmetry is declared, the function returns without perfor-
- ming any checks. Syntax:
- validate_sym(spin_system,bas)
- spin_system -Spinach spin system description object, with the

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ismember()`, `isfield()`, `any()`, `strcmp()`, `spins()`, `perm_group()`, `perm()`, `num2str()`, `get_coupling()`, `report()`, `iscell()`, `cellfun()`.
