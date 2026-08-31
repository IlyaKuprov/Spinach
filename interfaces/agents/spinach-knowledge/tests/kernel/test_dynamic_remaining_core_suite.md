# tests/kernel/test_dynamic_remaining_core_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_remaining_core_suite.m`
- Signature: `result=test_dynamic_remaining_core_suite()`
- Total lines: 287

## Purpose

Tests remaining deterministic utility helpers. Syntax: result=test_dynamic_remaining_core_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file also defines local helper function(s): `local_liouvillian_system()`, `local_dipolar_system()`, `local_isoswap_inputs()`, `gcp()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Remaining deterministic utility helpers\n')`.
- Lines 20-23: State the utility target of the test; implemented by `result=new_test_result('kernel/dynamic_remaining_core_suite', 'Remaining deterministic utility helpers', 'Small hand-written utility calls must preserve documented algeb…`.
- Lines 25-26: Build a minimal Liouville-space descriptor for algebraic utilities; implemented by `spin_system=local_liouvillian_system(3)`.
- Lines 28-29: Check adiabatic elimination against the exact Schur complement term; implemented by `L=[1 2 3;4 5 6;7 8 10]`.
- Lines 37-38: Check a textbook Clebsch-Gordan coefficient; implemented by `cg_triplet=cg_fast(1,0,1/2,1/2,1/2,-1/2)`.
- Lines 43-44: Check direct console and banner reporting into a file handle; implemented by `log_file=[tempname '.txt']`.
- Lines 57-58: Check polyadic text diagnostics with an explicit label; implemented by `polinfo(polyadic({{speye(2),sparse(1)}}),1,'labelled')`.
- Lines 61-62: Check summary traverses coordinate metadata through the silent report path; implemented by `spin_system.comp.nspins=1`.
- Lines 70-71: Check impound cell packaging preserves values and types; implemented by `payload=impound(17,'spinach',{speye(2)})`.
- Lines 77-78: Check dipolar coupling for two spins one Angstrom apart on the X axis; implemented by `dip_system=local_dipolar_system()`.
- Lines 89-90: Check isotope swapping rescales spin-pair couplings and wipes quadratic terms; implemented by `[sys,inter]=local_isoswap_inputs()`.
- Lines 102-103: Check interaction-representation order zero against the perturbation block; implemented by `spin_system=local_liouvillian_system(2)`.
- Lines 110-111: Check finite-difference Hessian construction on a constant 3D field; implemented by `H=fdhess(ones(5,6,7),3)`.
- Lines 122-123: Check sinkhole column removal in spherical-tensor Liouville formalism; implemented by `spin_system=local_liouvillian_system(3)`.
- Lines 129-130: Check the normalised Lorentzian branch analytically; implemented by `x_axis=[-1 0 1]`.
- Lines 136-137: Check the narrow-linewidth Lorentzian peak at exact centre; implemented by `tiny_fwhm=1e-308`.
- Lines 144-145: Check the narrow-linewidth MEX segment branch at exact offsets; implemented by `if exist('lorentzcon','file')==3`.
- Lines 152-153: Check the normalised Gaussian branch analytically; implemented by `x_axis=[-1 0 1]`.

### Control flow inferred from the code

- Line 114: `for` loop over `n=1:3`.
- Line 115: `for` loop over `k=1:3`.
- Line 145: conditional branch on `exist('lorentzcon','file')==3`.
- Line 170: `for` loop over `n=1:numel(x_axis)`.
- Line 212: `for` loop over `n=1:t2.nsteps`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/dynamic_remaining_core_suite', 'Remaining deterministic utility helpers', 'Small hand-written utility calls must preserve documented algeb…`.
- Lines 26: computes `spin_system` using `spin_system=local_liouvillian_system(3)`.
- Lines 29: computes `L` using `L=[1 2 3;4 5 6;7 8 10]`.
- Lines 30: computes `[L_slow,R_extra]` using `[L_slow,R_extra]=adelim(spin_system,L,3,[1 2])`.
- Lines 31: computes `R_ref` using `R_ref=1i*[3;6]*(1/10)*[7 8]`.
- Lines 38: computes `cg_triplet` using `cg_triplet=cg_fast(1,0,1/2,1/2,1/2,-1/2)`.
- Lines 40-41: computes `1e-12,1e-12, 'two spin-half states couple to the triplet M` using `1e-12,1e-12, 'two spin-half states couple to the triplet M=0 state with coefficient 1/sqrt(2)')`.
- Lines 41: computes `'two spin-half states couple to the triplet M` using `'two spin-half states couple to the triplet M=0 state with coefficient 1/sqrt(2)')`.
- Lines 44: computes `log_file` using `log_file=[tempname '.txt']`.
- Lines 45: computes `file_id` using `file_id=fopen(log_file,'w')`.
- Lines 46: computes `spin_system.sys.output` using `spin_system.sys.output=file_id`.
- Lines 50: computes `log_text` using `log_text=fileread(log_file)`.
- Lines 59: computes `result.messages{end+1}` using `result.messages{end+1}='PASS: polinfo labelled output -- direct polyadic diagnostic call completed'`.
- Lines 62: computes `spin_system.comp.nspins` using `spin_system.comp.nspins=1`.
- Lines 63: computes `spin_system.comp.isotopes` using `spin_system.comp.isotopes={'1H'}`.
- Lines 64: computes `spin_system.comp.labels` using `spin_system.comp.labels={'proton'}`.
- Lines 65: computes `spin_system.comp.mults` using `spin_system.comp.mults=2`.
- Lines 66: computes `spin_system.inter.coordinates` using `spin_system.inter.coordinates={[0 0 0]}`.

### Local helper functions

- Line 231: `local_liouvillian_system()` — `function spin_system=local_liouvillian_system(dim)`. Create a quiet spherical-tensor Liouville descriptor
  - Representative operation: `spin_system.sys.output='hush'`.
  - Representative operation: `spin_system.sys.enable={}`.
- Line 247: `local_dipolar_system()` — `function spin_system=local_dipolar_system()`. Create a quiet two-spin descriptor with one close internuclear vector
  - Representative operation: `spin_system.sys.output='hush'`.
  - Representative operation: `spin_system.sys.enable={}`.
- Line 266: `local_isoswap_inputs()` — `function [sys,inter]=local_isoswap_inputs()`. Create two-spin interaction structures with transferable and wiped terms
  - Representative operation: `sys.isotopes={'1H','13C'}`.
  - Representative operation: `inter.coupling.matrix=cell(2,2)`.
- Line 279: `local_ensure_pool()` — `function local_ensure_pool()`. Start a one-worker process pool for compact parfor utilities
  - Representative operation: `current_pool=gcp('nocreate')`.
  - Representative operation: `if isempty(current_pool)`.

## Outputs

- result -regression test result with explanatory messages
- The test checks block eliminations, text reporting, spin metadata,
- analytical line shapes, pumping terms, kite pruning, trajectory
- stitching, and small random-rotation diagnostics.

## Implementation structure

- Tests remaining deterministic utility helpers. Syntax:
- result=test_dynamic_remaining_core_suite()
- result -regression test result with explanatory messages
- The test checks block eliminations, text reporting, spin metadata,
- analytical line shapes, pumping terms, kite pruning, trajectory
- stitching, and small random-rotation diagnostics.
- Announce the test target
- State the utility target of the test
- Build a minimal Liouville-space descriptor for algebraic utilities
- Check adiabatic elimination against the exact Schur complement term
- Check a textbook Clebsch-Gordan coefficient
- Check direct console and banner reporting into a file handle

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `local_liouvillian_system()`, `adelim()`, `test_close()`, `cg_fast()`, `fopen()`, `report()`, `banner()`, `fclose()`, `fileread()`, `delete()`, `test_true()`, `contains()`, `polinfo()`, `polyadic()`, `speye()`.
