# kernel/utilities/kill_spin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/kill_spin.m`
- Signature: `spin_system=kill_spin(spin_system,hit_list)`
- Total lines: 194

## Purpose

Removes the specified spins from the spin_system structure and updates it accordingly. Syntax: spin_system=kill_spin(spin_system,hit_list)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 31-32: Check consistency; implemented by `grumble(spin_system,hit_list)`.
- Lines 34-35: Catch logical indexing; implemented by `if islogical(hit_list), hit_list=find(hit_list); end`.
- Lines 37-39: Inform the user; implemented by `report(spin_system,['removing ' num2str(numel(hit_list)) ' spins from the system '])`.
- Lines 41-42: Update isotope and particle type lists; implemented by `spin_system.comp.isotopes(hit_list)=[]`.
- Lines 45-46: Update the isotope list hash; implemented by `if isfield(spin_system.comp,'iso_hash')`.
- Lines 50-51: Update spin numbers; implemented by `spin_system.comp.nspins=spin_system.comp.nspins-numel(hit_list)`.
- Lines 53-54: Update labels list; implemented by `spin_system.comp.labels(hit_list)=[]`.
- Lines 56-57: Update multiplicities and magnetogyric ratios; implemented by `spin_system.comp.mults(hit_list)=[]`.
- Lines 60-61: Update base frequencies; implemented by `spin_system.inter.basefrqs(hit_list)=[]`.
- Lines 63-64: Update Zeeman tensor array; implemented by `spin_system.inter.zeeman.matrix(hit_list)=[]`.
- Lines 66-67: Update DD scaling multipliers; implemented by `spin_system.inter.zeeman.ddscal(hit_list)=[]`.
- Lines 69-70: Update giant spin Hamiltonian terms; implemented by `spin_system.inter.giant.coeff(hit_list)=[]`.
- Lines 72-73: Update coupling tensor array; implemented by `spin_system.inter.coupling.matrix(hit_list,:)=[]`.
- Lines 76-77: Update coordinates; implemented by `spin_system.inter.coordinates(hit_list)=[]`.
- Lines 79-80: Update proximity matrix; implemented by `spin_system.inter.proxmatrix(hit_list,:)=[]`.
- Lines 83-84: Update relaxation parameters; implemented by `if ~isempty(spin_system.rlx.r1_rates)`.
- Lines 109-110: Update scalar relaxation source spins; implemented by `if ~isempty(spin_system.rlx.srsk_sources)`.
- Lines 116-117: Update kinetics parameters; implemented by `for n=1:numel(spin_system.chem.parts)`.

### Control flow inferred from the code

- Line 35: conditional branch on `islogical(hit_list), hit_list=find(hit_list); end`.
- Line 46: conditional branch on `isfield(spin_system.comp,'iso_hash')`.
- Line 84: conditional branch on `~isempty(spin_system.rlx.r1_rates)`.
- Line 87: conditional branch on `~isempty(spin_system.rlx.r2_rates)`.
- Line 90: conditional branch on `~isempty(spin_system.rlx.lind_r1_rates)`.
- Line 93: conditional branch on `~isempty(spin_system.rlx.lind_r2_rates)`.
- Line 96: conditional branch on `~isempty(spin_system.rlx.srfk_mdepth)`.
- Line 100: conditional branch on `~isempty(spin_system.rlx.weiz_r1d)`.
- Line 104: conditional branch on `~isempty(spin_system.rlx.weiz_r2d)`.
- Line 110: conditional branch on `~isempty(spin_system.rlx.srsk_sources)`.
- Line 117: `for` loop over `n=1:numel(spin_system.chem.parts)`.
- Line 123: conditional branch on `~isempty(spin_system.chem.flux_rate)`.
- Line 132: conditional branch on `(~isempty(spin_system.chem.rp_rates))&&(numel(spin_system.chem.rp_electrons)<2)`.
- Line 137: conditional branch on `isfield(spin_system,'bas')`.

### Key state/data transformations

- Lines 42: computes `spin_system.comp.isotopes(hit_list)` using `spin_system.comp.isotopes(hit_list)=[]`.
- Lines 43: computes `spin_system.comp.types(hit_list)` using `spin_system.comp.types(hit_list)=[]`.
- Lines 47: computes `spin_system.comp.iso_hash` using `spin_system.comp.iso_hash=md5_hash(spin_system.comp.isotopes)`.
- Lines 51: computes `spin_system.comp.nspins` using `spin_system.comp.nspins=spin_system.comp.nspins-numel(hit_list)`.
- Lines 54: computes `spin_system.comp.labels(hit_list)` using `spin_system.comp.labels(hit_list)=[]`.
- Lines 57: computes `spin_system.comp.mults(hit_list)` using `spin_system.comp.mults(hit_list)=[]`.
- Lines 58: computes `spin_system.inter.gammas(hit_list)` using `spin_system.inter.gammas(hit_list)=[]`.
- Lines 61: computes `spin_system.inter.basefrqs(hit_list)` using `spin_system.inter.basefrqs(hit_list)=[]`.
- Lines 64: computes `spin_system.inter.zeeman.matrix(hit_list)` using `spin_system.inter.zeeman.matrix(hit_list)=[]`.
- Lines 67: computes `spin_system.inter.zeeman.ddscal(hit_list)` using `spin_system.inter.zeeman.ddscal(hit_list)=[]`.
- Lines 70: computes `spin_system.inter.giant.coeff(hit_list)` using `spin_system.inter.giant.coeff(hit_list)=[]`.
- Lines 73: computes `spin_system.inter.coupling.matrix(hit_list,:)` using `spin_system.inter.coupling.matrix(hit_list,:)=[]`.
- Lines 74: computes `spin_system.inter.coupling.matrix(:,hit_list)` using `spin_system.inter.coupling.matrix(:,hit_list)=[]`.
- Lines 77: computes `spin_system.inter.coordinates(hit_list)` using `spin_system.inter.coordinates(hit_list)=[]`.
- Lines 80: computes `spin_system.inter.proxmatrix(hit_list,:)` using `spin_system.inter.proxmatrix(hit_list,:)=[]`.
- Lines 81: computes `spin_system.inter.proxmatrix(:,hit_list)` using `spin_system.inter.proxmatrix(:,hit_list)=[]`.
- Lines 85: computes `spin_system.rlx.r1_rates(hit_list)` using `spin_system.rlx.r1_rates(hit_list)=[]`.
- Lines 88: computes `spin_system.rlx.r2_rates(hit_list)` using `spin_system.rlx.r2_rates(hit_list)=[]`.

### Local helper functions

- Line 173: `grumble()` — `function grumble(spin_system,hit_list)`.
  - Representative operation: `if islogical(hit_list)`.
  - Representative operation: `if(numel(hit_list)~=spin_system.comp.nspins)`.

## Parameters / inputs

- spin_system -primary Spinach data structure
- hit_list -a vector of integers or a logical
- vector giving the numbers of spins
- to be removed from the system

## Outputs

- spin_system -the data structure with the indica-
- ted spins and all dependent infor-
- mation (basis, assumptions) removed
- Notes: basis, connectivity, symmetry, and assumption information
- is destroyed by this function; you would need to call the
- basis.m and assume.m functions again.

## Implementation structure

- Removes the specified spins from the spin_system structure and
- updates it accordingly. Syntax:
- spin_system=kill_spin(spin_system,hit_list)
- spin_system -primary Spinach data structure
- hit_list -a vector of integers or a logical
- vector giving the numbers of spins
- to be removed from the system
- spin_system -the data structure with the indica-
- ted spins and all dependent infor-
- mation (basis, assumptions) removed
- is destroyed by this function; you would need to call the
- basis.m and assume.m functions again.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `islogical()`, `report()`, `num2str()`, `isfield()`, `md5_hash()`, `srsk_spins()`, `false()`, `subsystem_idx()`, `true()`, `reacting_spins()`, `rmfield()`, `any()`.
