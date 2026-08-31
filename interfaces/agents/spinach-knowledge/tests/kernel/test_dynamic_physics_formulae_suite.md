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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Physical formula utilities\n')`.
- Lines 20-23: State the physical formula target of the test; implemented by `result=new_test_result('kernel/dynamic_physics_formulae_suite', 'Physical formula utilities', 'closed-form physical helper formulae must match their analytic definitions…`.
- Lines 25-26: Check Clebsch-Gordan reduction of two spin-half irreps; implemented by `[mults,proj]=add_spins(1/2,1/2)`.
- Lines 35-36: Check point-dipole coupling for a one-Angstrom z-axis displacement; implemented by `hbar=1.054571628e-34`.
- Lines 47-48: Check point electron-nucleus hyperfine tensor for a z-axis displacement; implemented by `A=xyz2hfc([0 0 0],[0 0 1],'1H')`.
- Lines 53-54: Check exponential drop values at exact quartering points; implemented by `rate=log(4)`.
- Lines 59-60: Check skew-normal density reduces to the normal density at zero skew; implemented by `x=[-1 0 1]`.
- Lines 66-67: Check oscillator coordinate grid and coordinate operator; implemented by `parameters.frc_cnst=4`.
- Lines 80-81: Check one-dimensional hydrodynamic derivative construction; implemented by `hydro_spin_system.sys.enable={'polyadic'}`.
- Lines 92-93: Check spherical-tensor to Zeeman projection metadata on one spin-half; implemented by `spin_system.comp.mults=2`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/dynamic_physics_formulae_suite', 'Physical formula utilities', 'closed-form physical helper formulae must match their analytic definitions…`.
- Lines 26: computes `[mults,proj]` using `[mults,proj]=add_spins(1/2,1/2)`.
- Lines 27: computes `P` using `P=[proj{1} proj{2}]`.
- Lines 36: computes `hbar` using `hbar=1.054571628e-34`.
- Lines 37: computes `mu0` using `mu0=4*pi*1e-7`.
- Lines 38: computes `[d,alp,bet,gam,M]` using `[d,alp,bet,gam,M]=xyz2dd([0 0 0],[0 0 1],'1H','13C')`.
- Lines 39: computes `d_ref` using `d_ref=spin('1H')*spin('13C')*hbar*mu0/(4*pi*(1e-10)^3)`.
- Lines 48: computes `A` using `A=xyz2hfc([0 0 0],[0 0 1],'1H')`.
- Lines 49: computes `C` using `C=1e4*spin('1H')*hbar*mu0/(4*pi*(1e-10)^3)`.
- Lines 54: computes `rate` using `rate=log(4)`.
- Lines 55: computes `drop` using `drop=expdrop(5,1,2,3,rate)`.
- Lines 60: computes `x` using `x=[-1 0 1]`.
- Lines 61: computes `p` using `p=snormpdf(x,0,1,0)`.
- Lines 62: computes `p_ref` using `p_ref=exp(-(x.^2)/2)/sqrt(2*pi)`.
- Lines 67: computes `parameters.frc_cnst` using `parameters.frc_cnst=4`.
- Lines 68: computes `parameters.par_mass` using `parameters.par_mass=2`.
- Lines 69: computes `parameters.grv_cnst` using `parameters.grv_cnst=0`.
- Lines 70: computes `parameters.n_points` using `parameters.n_points=11`.

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
