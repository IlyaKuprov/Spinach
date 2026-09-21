# tests/kernel/test_dynamic_physics_formulae_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_physics_formulae_suite.m`
- Signature: `result=test_dynamic_physics_formulae_suite()`
- Total lines: 102

## Purpose

Tests deterministic physical formula utility helpers. Syntax: result=test_dynamic_physics_formulae_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.

## Outputs

- result -regression test result with explanatory messages
- The test checks spin addition, point-dipole tensors, hyperfine tensors,
- exponential drops, skew-normal densities, oscillator grids, hydrodynamic
- derivative construction, and spherical-tensor projection metadata.

## Implementation structure

- Tests deterministic physical formula utility helpers. Syntax:
- result=test_dynamic_physics_formulae_suite()
- result -regression test result with explanatory messages
- The test checks spin addition, point-dipole tensors, hyperfine tensors,
- exponential drops, skew-normal densities, oscillator grids, hydrodynamic
- derivative construction, and spherical-tensor projection metadata.
- Announce the test target
- State the physical formula target of the test
- Check Clebsch-Gordan reduction of two spin-half irreps
- Check point-dipole coupling for a one-Angstrom z-axis displacement
- Check point electron-nucleus hyperfine tensor for a z-axis displacement
- Check exponential drop values at exact quartering points

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `add_spins()`, `test_true()`, `isequal()`, `test_close()`, `xyz2dd()`, `spin()`, `xyz2hfc()`, `expdrop()`, `snormpdf()`, `oscillator()`, `hydrodynamics()`, `fdmat()`, `inflate()`, `sphten2zeeman()`.
