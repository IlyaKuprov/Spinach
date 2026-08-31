# tests/kernel/test_dynamic_chem_geometry_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_chem_geometry_suite.m`
- Signature: `result=test_dynamic_chem_geometry_suite()`
- Total lines: 130

## Purpose

Tests deterministic chemistry and geometry utility helpers. Syntax: result=test_dynamic_chem_geometry_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file also defines local helper function(s): `local_geometry_spin_system()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Chemistry and geometry utilities\n')`.
- Lines 19-22: State the chemistry and geometry target of the test; implemented by `result=new_test_result('kernel/dynamic_chem_geometry_suite', 'Chemistry and geometry utilities', 'small chemistry and geometry helpers must preserve exact coordinate, te…`.
- Lines 24-25: Build a small spin-system descriptor for metadata helpers; implemented by `spin_system=local_geometry_spin_system()`.
- Lines 27-28: Check simple cubic lattice construction and periodic vectors; implemented by `[sys,inter]=cubic_lattice('13C',2,2)`.
- Lines 36-37: Check a signed right-angle dihedral from four Cartesian points; implemented by `phi=dihedral([1 0 0],[0 0 0],[0 1 0],[0 1 1])`.
- Lines 41-42: Check Cartesian point-cloud binning on a regular grid; implemented by `coords=[0.25 0.25 0.25;0.75 0.25 0.25;1.5 0.5 0.5]`.
- Lines 48-49: Check nearest-neighbour lookup by Cartesian coordinate distance; implemented by `[spin_idx,dist]=nearest_spin(spin_system,1)`.
- Lines 55-57: Check chemical-substance ownership lookup; implemented by `result=test_close(result,'which_subst first part',which_subst(spin_system,[1 2]),1,1e-15,1e-15, 'spins one and two both belong to the first chemical part')`.
- Lines 61-62: Check coupling extraction from forward and backward tensor entries; implemented by `coupling=get_coupling(spin_system,1,2)`.
- Lines 66-67: Check chemical shifts and Zeeman offsets from isotropic tensor parts; implemented by `[cs_ppm,cs_hz]=chemshifts(spin_system)`.
- Lines 75-76: Check g-tensor extraction from the stored Zeeman scaling tensor; implemented by `spin_system.inter.zeeman.ddscal{1}=2*eye(3)`.
- Lines 82-83: Check replacement of isotropic tensor components without changing anisotropy; implemented by `input_tensors={diag([1 2 6]),eye(3)}`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/dynamic_chem_geometry_suite', 'Chemistry and geometry utilities', 'small chemistry and geometry helpers must preserve exact coordinate, te…`.
- Lines 25: computes `spin_system` using `spin_system=local_geometry_spin_system()`.
- Lines 28: computes `[sys,inter]` using `[sys,inter]=cubic_lattice('13C',2,2)`.
- Lines 37: computes `phi` using `phi=dihedral([1 0 0],[0 0 0],[0 1 0],[0 1 1])`.
- Lines 42: computes `coords` using `coords=[0.25 0.25 0.25;0.75 0.25 0.25;1.5 0.5 0.5]`.
- Lines 43: computes `density` using `density=xyz2pd(coords,[0 1],[0 1],[0 1],2,2,2)`.
- Lines 44: computes `density_ref` using `density_ref=zeros(2,2,2); density_ref(1,1,1)=1; density_ref(2,1,1)=1`.
- Lines 49: computes `[spin_idx,dist]` using `[spin_idx,dist]=nearest_spin(spin_system,1)`.
- Lines 62: computes `coupling` using `coupling=get_coupling(spin_system,1,2)`.
- Lines 67: computes `[cs_ppm,cs_hz]` using `[cs_ppm,cs_hz]=chemshifts(spin_system)`.
- Lines 76: computes `spin_system.inter.zeeman.ddscal{1}` using `spin_system.inter.zeeman.ddscal{1}=2*eye(3)`.
- Lines 77-78: computes `g_ref` using `g_ref=-spin_system.inter.zeeman.ddscal{1}*spin_system.inter.gammas(1)* spin_system.tols.hbar/spin_system.tols.muB`.
- Lines 83: computes `input_tensors` using `input_tensors={diag([1 2 6]),eye(3)}`.
- Lines 84: computes `shifted` using `shifted=shift_iso(input_tensors,1,10)`.

### Local helper functions

- Line 93: `local_geometry_spin_system()` — `function spin_system=local_geometry_spin_system()`. Create quiet system output settings
  - Representative operation: `spin_system.sys.output='hush'`.
  - Representative operation: `spin_system.sys.disable={}`.

## Outputs

- result -regression test result with explanatory messages
- The test checks lattice construction, geometric measurements, coupling
- extraction, chemical shifts, nearest-neighbour lookup, and tensor helpers.

## Implementation structure

- Tests deterministic chemistry and geometry utility helpers. Syntax:
- result=test_dynamic_chem_geometry_suite()
- result -regression test result with explanatory messages
- The test checks lattice construction, geometric measurements, coupling
- extraction, chemical shifts, nearest-neighbour lookup, and tensor helpers.
- Announce the test target
- State the chemistry and geometry target of the test
- Build a small spin-system descriptor for metadata helpers
- Check simple cubic lattice construction and periodic vectors
- Check a signed right-angle dihedral from four Cartesian points
- Check Cartesian point-cloud binning on a regular grid
- Check nearest-neighbour lookup by Cartesian coordinate distance

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `local_geometry_spin_system()`, `cubic_lattice()`, `test_true()`, `all()`, `strcmp()`, `test_close()`, `isequal()`, `dihedral()`, `xyz2pd()`, `density_ref()`, `nearest_spin()`, `which_subst()`, `get_coupling()`, `chemshifts()`, `cs_ppm()`.
