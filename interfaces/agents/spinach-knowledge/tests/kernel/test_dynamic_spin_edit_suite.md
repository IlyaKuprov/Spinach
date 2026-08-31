# tests/kernel/test_dynamic_spin_edit_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_spin_edit_suite.m`
- Signature: `result=test_dynamic_spin_edit_suite()`
- Total lines: 270

## Purpose

Tests deterministic spin-system editing support utilities. Syntax: result=test_dynamic_spin_edit_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The file also defines local helper function(s): `local_merge_fails()`, `local_edit_spin_system()`, `local_assumed_spin_system()`, `local_merge_parts()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Spin-system editing utilities\n')`.
- Lines 19-22: State the spin editing target of the test; implemented by `result=new_test_result('kernel/dynamic_spin_edit_suite', 'Spin-system editing utilities', 'local spin-system editing helpers must update dependent metadata without touch…`.
- Lines 24-25: Build a small but structurally complete spin-system descriptor; implemented by `spin_system=local_edit_spin_system()`.
- Lines 27-28: Check numeric spin removal updates identities, labels, coordinates, and parts; implemented by `trimmed=kill_spin(spin_system,2)`.
- Lines 44-45: Check chemical part renumbering across multiple subsystems; implemented by `multipart=spin_system; multipart.chem.parts={[1 2],3}`.
- Lines 50-51: Check destruction of stale basis, connectivity, symmetry, and assumption data; implemented by `stale=spin_system; stale.bas.formalism='sphten-liouv'`.
- Lines 65-66: Check logical spin removal follows the same path; implemented by `logical_trimmed=kill_spin(spin_system,[false true false])`.
- Lines 71-72: Check dilute-isotope subsystem generation through kill_spin; implemented by `subsystems=dilute(spin_system,'13C',1)`.
- Lines 78-79: Check assumption overrides by numeric and isotope specifications; implemented by `assumed=local_assumed_spin_system()`.
- Lines 89-90: Check merging of independent Spinach input structure fragments; implemented by `[sys_parts,inter_parts]=local_merge_parts()`.
- Lines 117-118: Check column-oriented subsystem lists and partless chemistry; implemented by `[sys_parts,inter_parts]=local_merge_parts()`.
- Lines 132-133: Check that non-extensive differences and malformed inputs are refused; implemented by `[sys_parts,inter_parts]=local_merge_parts()`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/dynamic_spin_edit_suite', 'Spin-system editing utilities', 'local spin-system editing helpers must update dependent metadata without touch…`.
- Lines 25: computes `spin_system` using `spin_system=local_edit_spin_system()`.
- Lines 28: computes `trimmed` using `trimmed=kill_spin(spin_system,2)`.
- Lines 45: computes `multipart` using `multipart=spin_system; multipart.chem.parts={[1 2],3}`.
- Lines 51: computes `stale` using `stale=spin_system; stale.bas.formalism='sphten-liouv'`.
- Lines 52: computes `stale.inter.conmatrix` using `stale.inter.conmatrix=logical(speye(3))`.
- Lines 53: computes `stale.comp.sym_group` using `stale.comp.sym_group={'S2'}; stale.comp.sym_spins={[2 3]}; stale.comp.sym_a1g_only=true()`.
- Lines 54: computes `stale.inter.assumptions` using `stale.inter.assumptions='nmr'`.
- Lines 55: computes `stale.inter.zeeman.strength` using `stale.inter.zeeman.strength={'secular','secular','secular'}`.
- Lines 56: computes `stale.inter.giant.strength` using `stale.inter.giant.strength={[],[],[]}`.
- Lines 57: computes `stale.inter.coupling.strength` using `stale.inter.coupling.strength=cell(3,3)`.
- Lines 66: computes `logical_trimmed` using `logical_trimmed=kill_spin(spin_system,[false true false])`.
- Lines 72: computes `subsystems` using `subsystems=dilute(spin_system,'13C',1)`.
- Lines 73: computes `carbon_counts` using `carbon_counts=cellfun(@(x)nnz(strcmp(x.comp.isotopes,'13C')),subsystems)`.
- Lines 74: computes `spin_counts` using `spin_counts=cellfun(@(x)x.comp.nspins,subsystems)`.
- Lines 79: computes `assumed` using `assumed=local_assumed_spin_system()`.
- Lines 90: computes `[sys_parts,inter_parts]` using `[sys_parts,inter_parts]=local_merge_parts()`.
- Lines 91: computes `[sys,inter]` using `[sys,inter]=merge_inp(sys_parts,inter_parts)`.

### Local helper functions

- Line 158: `local_merge_fails()` — `function failed=local_merge_fails(sys_parts,inter_parts)`.
  - Representative operation: `try`.
  - Representative operation: `merge_inp(sys_parts,inter_parts); failed=false()`.
- Line 167: `local_edit_spin_system()` — `function spin_system=local_edit_spin_system()`. Create quiet system output settings
  - Representative operation: `spin_system.sys.output='hush'`.
  - Representative operation: `spin_system.sys.enable={}`.
- Line 211: `local_assumed_spin_system()` — `function spin_system=local_assumed_spin_system()`. Start from the edit fixture and add assumption fields
  - Representative operation: `spin_system=local_edit_spin_system()`.
  - Representative operation: `spin_system.inter.zeeman.strength={'secular','secular','secular'}`.
- Line 222: `local_merge_parts()` — `function [sys_parts,inter_parts]=local_merge_parts()`. Build first subsystem input structures
  - Representative operation: `sys_parts{1}.magnet=14.1`.
  - Representative operation: `sys_parts{1}.isotopes={'1H'}`.

## Outputs

- result -regression test result with explanatory messages
- The test checks spin removal, dilute-isotope subsystem generation,
- assumption overrides, and merging of Spinach input structures.

## Implementation structure

- Tests deterministic spin-system editing support utilities. Syntax:
- result=test_dynamic_spin_edit_suite()
- result -regression test result with explanatory messages
- The test checks spin removal, dilute-isotope subsystem generation,
- assumption overrides, and merging of Spinach input structures.
- Announce the test target
- State the spin editing target of the test
- Build a small but structurally complete spin-system descriptor
- Check numeric spin removal updates identities, labels, coordinates, and parts
- Check chemical part renumbering across multiple subsystems
- Check destruction of stale basis, connectivity, symmetry, and assumption data
- Check logical spin removal follows the same path

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `local_edit_spin_system()`, `kill_spin()`, `test_true()`, `isequal()`, `strcmp()`, `md5_hash()`, `logical()`, `speye()`, `true()`, `isfield()`, `dilute()`, `cellfun()`, `nnz()`, `all()`, `local_assumed_spin_system()`.
