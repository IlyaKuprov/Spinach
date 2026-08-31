# kernel/utilities/merge_inp.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/merge_inp.m`
- Signature: `[sys,inter]=merge_inp(sys_parts,inter_parts)`
- Total lines: 443

## Purpose

Merges multiple sys and inter structures into one. Useful for setting up chemical kinetics simulations where the molecules come from different DFT calculations. Syntax: [sys,inter]=merge_inp(sys_parts,inter_parts)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `group_check()`, `group_gate()`, `strip()`, `merge_like_magnet()`, `merge_like_isotopes()`, `merge_like_coords()`, `merge_like_rates()`, `merge_like_couplings()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 40-41: Check consistency; implemented by `grumble(sys_parts,inter_parts)`.
- Lines 43-44: Count the spins in each subsystem; implemented by `part_sizes=cellfun(@(x)numel(x.isotopes),sys_parts)`.
- Lines 46-47: Create structure stubs; implemented by `sys.stub=1; inter.stub=1`.
- Lines 49-51: Fields containing common values that must be equal; implemented by `sys_common={'magnet','output','scratch','enable', 'disable','tols','parallel','parprops'}`.
- Lines 56-57: Fields containing row cell arrays; implemented by `[sys,sys_parts]=merge_like_isotopes(sys,sys_parts,'isotopes')`.
- Lines 60-65: Fields containing common values that must be equal; implemented by `inter_common={'relaxation','rlx_keep','rlx_dfs','equilibrium', 'temperature','damp_rate','srfk_tau_c', 'nz_shift','nz_onshell', 'weiz_r1e','weiz_r2e','weiz_r1n','weiz_r2…`.
- Lines 70-72: Fields that are per chemical species when a species split is present; implemented by `if group_check(inter_parts,'chem')&& all(cellfun(@(x)isfield(x.chem,'parts'),inter_parts))`.
- Lines 80-81: Fields containing column cell arrays; implemented by `[inter,inter_parts]=merge_like_coords(inter,inter_parts,'coordinates')`.
- Lines 83-84: Fields containing per-spin relaxation rates; implemented by `inter_rates={'r1_rates','r2_rates','lind_r1_rates','lind_r2_rates'}`.
- Lines 89-90: Fields containing square arrays over spin pairs; implemented by `inter_square={'srfk_mdepth','weiz_r1d','weiz_r2d'}`.
- Lines 95-96: Fields containing spin index vectors; implemented by `[inter,inter_parts]=merge_like_sources(inter,inter_parts,'srsk_sources',part_sizes)`.
- Lines 98-99: Fields containing cell arrays of spin index vectors; implemented by `[inter,inter_parts]=merge_like_parts(inter,inter_parts,'ignore',part_sizes)`.
- Lines 101-102: Zeeman interaction specifications; implemented by `if group_check(inter_parts,'zeeman')`.
- Lines 114-115: Coupling specifications; implemented by `if group_check(inter_parts,'coupling')`.
- Lines 127-128: Giant spin interaction specifications; implemented by `if group_check(inter_parts,'giant')`.
- Lines 138-139: Magnetic susceptibility centre specifications; implemented by `if group_check(inter_parts,'suscept')`.
- Lines 149-150: Chemical process specifications; implemented by `if group_check(inter_parts,'chem')`.
- Lines 166-167: Catch unhandled subfields; implemented by `for n=1:numel(sys_parts)`.

### Control flow inferred from the code

- Line 52: `for` loop over `n=1:numel(sys_common)`.
- Line 66: `for` loop over `n=1:numel(inter_common)`.
- Line 71: conditional branch on `group_check(inter_parts,'chem')&&`.
- Line 85: `for` loop over `n=1:numel(inter_rates)`.
- Line 91: `for` loop over `n=1:numel(inter_square)`.
- Line 102: conditional branch on `group_check(inter_parts,'zeeman')`.
- Line 115: conditional branch on `group_check(inter_parts,'coupling')`.
- Line 128: conditional branch on `group_check(inter_parts,'giant')`.
- Line 139: conditional branch on `group_check(inter_parts,'suscept')`.
- Line 150: conditional branch on `group_check(inter_parts,'chem')`.
- Line 167: `for` loop over `n=1:numel(sys_parts)`.
- Line 168: conditional branch on `~isempty(fieldnames(sys_parts{n}))`.
- Line 172: `for` loop over `n=1:numel(inter_parts)`.
- Line 173: conditional branch on `~isempty(fieldnames(inter_parts{n}))`.

### Key state/data transformations

- Lines 44: computes `part_sizes` using `part_sizes=cellfun(@(x)numel(x.isotopes),sys_parts)`.
- Lines 47: computes `sys.stub` using `sys.stub=1; inter.stub=1`.
- Lines 50-51: computes `sys_common` using `sys_common={'magnet','output','scratch','enable', 'disable','tols','parallel','parprops'}`.
- Lines 53: computes `[sys,sys_parts]` using `[sys,sys_parts]=merge_like_magnet(sys,sys_parts,sys_common{n})`.
- Lines 61-65: computes `inter_common` using `inter_common={'relaxation','rlx_keep','rlx_dfs','equilibrium', 'temperature','damp_rate','srfk_tau_c', 'nz_shift','nz_onshell', 'weiz_r1e','weiz_r2e','weiz_r1n','weiz_r2…`.
- Lines 67: computes `[inter,inter_parts]` using `[inter,inter_parts]=merge_like_magnet(inter,inter_parts,inter_common{n})`.
- Lines 84: computes `inter_rates` using `inter_rates={'r1_rates','r2_rates','lind_r1_rates','lind_r2_rates'}`.
- Lines 90: computes `inter_square` using `inter_square={'srfk_mdepth','weiz_r1d','weiz_r2d'}`.
- Lines 103: computes `zeeman_parts` using `zeeman_parts=strip(inter_parts,'zeeman'); zeeman.stub=1`.
- Lines 104: computes `[zeeman,zeeman_parts]` using `[zeeman,zeeman_parts]=merge_like_isotopes(zeeman,zeeman_parts,'matrix')`.
- Lines 109: computes `zeeman` using `zeeman=rmfield(zeeman,'stub'); inter.zeeman=zeeman`.
- Lines 110-111: computes `inter_parts` using `inter_parts=cellfun(@(x)rmfield(x,'zeeman'), inter_parts,'UniformOutput',false)`.
- Lines 116: computes `coupling_parts` using `coupling_parts=strip(inter_parts,'coupling'); coupling.stub=1`.
- Lines 117: computes `[coupling,coupling_parts]` using `[coupling,coupling_parts]=merge_like_couplings(coupling,coupling_parts,'matrix')`.
- Lines 122: computes `coupling` using `coupling=rmfield(coupling,'stub'); inter.coupling=coupling`.
- Lines 129: computes `giant_parts` using `giant_parts=strip(inter_parts,'giant'); giant.stub=1`.
- Lines 130: computes `[giant,giant_parts]` using `[giant,giant_parts]=merge_like_isotopes(giant,giant_parts,'coeff')`.
- Lines 133: computes `giant` using `giant=rmfield(giant,'stub'); inter.giant=giant`.

### Local helper functions

- Line 184: `group_check()` — `function group_present=group_check(parts_array,group_name)`. Catch unhandled subfields inside a nested group
  - Representative operation: `present=cellfun(@(x)isfield(x,group_name),parts_array)`.
  - Representative operation: `if any(present)&&(~all(present))`.
- Line 193: `group_gate()` — `function group_gate(group_parts,group_name)`. Strip a common field from an array of structures
  - Representative operation: `for n=1:numel(group_parts)`.
  - Representative operation: `if ~isempty(fieldnames(group_parts{n}))`.
- Line 202: `strip()` — `function parts_array=strip(parts_array,field_name)`. Merge subsystem fields assuming they hold an optional common value
  - Representative operation: `parts_array=cellfun(@(x)x.(field_name),parts_array,'UniformOutput',false)`.
- Line 207: `merge_like_magnet()` — `function [spec,spec_parts]=merge_like_magnet(spec,spec_parts,field_name)`.
  - Representative operation: `if isfield(spec_parts{1},field_name)`.
  - Representative operation: `spec.(field_name)=spec_parts{1}.(field_name)`.
- Line 235: `merge_like_isotopes()` — `function [spec,spec_parts]=merge_like_isotopes(spec,spec_parts,field_name)`.
  - Representative operation: `if isfield(spec_parts{1},field_name)`.
  - Representative operation: `spec.(field_name)=spec_parts{1}.(field_name)`.
- Line 263: `merge_like_coords()` — `function [spec,spec_parts]=merge_like_coords(spec,spec_parts,field_name)`.
  - Representative operation: `if isfield(spec_parts{1},field_name)`.
  - Representative operation: `spec.(field_name)=spec_parts{1}.(field_name)`.
- Line 291: `merge_like_rates()` — `function [spec,spec_parts]=merge_like_rates(spec,spec_parts,field_name)`.
  - Representative operation: `if isfield(spec_parts{1},field_name)`.
  - Representative operation: `spec.(field_name)=spec_parts{1}.(field_name)(:)`.
- Line 319: `merge_like_couplings()` — `function [spec,spec_parts]=merge_like_couplings(spec,spec_parts,field_name)`.
  - Representative operation: `if isfield(spec_parts{1},field_name)`.
  - Representative operation: `spec.(field_name)=spec_parts{1}.(field_name)`.

## Parameters / inputs

- sys_parts -a cell array of sys structures
- to be merged
- inter_parts -a cell array of inter structures
- to be merged

## Outputs

- sys -resulting sys structure
- inter -resulting inter structure
- Note: extensive fields are concatenated with spin and chemical
- subsystem indices offset as appropriate; non-extensive
- fields (magnet, temperature, relaxation settings, etc.)
- must be identical in all subsystems that supply them, and
- any difference is treated as an error. Nested groups
- (zeeman, coupling, giant, suscept, chem) and every field
- must be present in all subsystems or in none. An error is
- thrown for unhandled subfields. Coordinates and suscepti-
- bility centres from all subsystems are assumed to refer
- to one common frame of reference; spin index lists are
- returned as row vectors.

## Implementation structure

- Merges multiple sys and inter structures into one. Useful for
- setting up chemical kinetics simulations where the molecules
- come from different DFT calculations. Syntax:
- [sys,inter]=merge_inp(sys_parts,inter_parts)
- sys_parts -a cell array of sys structures
- to be merged
- inter_parts -a cell array of inter structures
- sys -resulting sys structure
- inter -resulting inter structure
- Note: extensive fields are concatenated with spin and chemical
- subsystem indices offset as appropriate; non-extensive
- fields (magnet, temperature, relaxation settings, etc.)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `cellfun()`, `merge_like_magnet()`, `merge_like_isotopes()`, `group_check()`, `all()`, `isfield()`, `merge_like_coords()`, `merge_like_rates()`, `merge_like_couplings()`, `merge_like_sources()`, `merge_like_parts()`, `strip()`, `group_gate()`, `rmfield()`, `fieldnames()`.
