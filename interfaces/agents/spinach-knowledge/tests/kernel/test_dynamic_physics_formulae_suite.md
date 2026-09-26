# tests/kernel/test_dynamic_physics_formulae_suite.m

- Signature: `result=test_dynamic_physics_formulae_suite()`

## Purpose

Tests deterministic physical formula utility helpers. Syntax: result=test_dynamic_physics_formulae_suite()

## Physical / mathematical content

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
