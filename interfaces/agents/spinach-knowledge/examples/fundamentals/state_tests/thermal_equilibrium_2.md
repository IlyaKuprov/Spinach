# examples/fundamentals/state_tests/thermal_equilibrium_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/state_tests/thermal_equilibrium_2.m`
- Signature: `thermal_equilibrium_2()`
- Total lines: 47

## Purpose

Thermal equilibrium states, using all the different formalisms supported by Spinach kernel.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 8-9: Spin system parameters; implemented by `sys.magnet=14.1`.
- Lines 20-21: Get thermal equilibrium states at finite temperature; implemented by `formalisms={'zeeman-hilb','zeeman-liouv','sphten-liouv'}`.
- Lines 31-32: Zeeman-liouv to zeeman-hilb; implemented by `rho{2}=reshape(rho{2},[96 96])`.
- Lines 34-35: Sphten-liouv to zeeman-hilb; implemented by `mult_factor=prod(spin_system.comp.mults)`.
- Lines 39-40: Check the differences; implemented by `if (norm(rho{1}-rho{2},2)>1e-8)||(norm(rho{2}-rho{3},2)>1e-8)`.

### Control flow inferred from the code

- Line 23: `for` loop over `n=1:numel(formalisms)`.
- Line 40: conditional branch on `(norm(rho{1}-rho{2},2)>1e-8)||(norm(rho{2}-rho{3},2)>1e-8)`.

### Key state/data transformations

- Lines 9: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 10: computes `sys.isotopes` using `sys.isotopes={'E8','1H','14N','15N'}`.
- Lines 11: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={2.002319 1.0 2.0 3.0}`.
- Lines 12: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(length(sys.isotopes))`.
- Lines 13: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=1e6`.
- Lines 14: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=1e6`.
- Lines 15: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=1e3`.
- Lines 16: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}=1e3`.
- Lines 17: computes `inter.coupling.scalar{1,4}` using `inter.coupling.scalar{1,4}=1e2`.
- Lines 18: computes `inter.temperature` using `inter.temperature=4.2`.
- Lines 21: computes `formalisms` using `formalisms={'zeeman-hilb','zeeman-liouv','sphten-liouv'}`.
- Lines 22: computes `rho` using `rho=cell(1,numel(formalisms))`.
- Lines 24: computes `bas.formalism` using `bas.formalism=formalisms{n}`.
- Lines 25: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 26: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 28: computes `rho{n}` using `rho{n}=equilibrium(spin_system)`.
- Lines 32: computes `rho{2}` using `rho{2}=reshape(rho{2},[96 96])`.
- Lines 35: computes `mult_factor` using `mult_factor=prod(spin_system.comp.mults)`.

## Implementation structure

- Thermal equilibrium states, using all the different formalisms
- supported by Spinach kernel.
- Spin system parameters
- Get thermal equilibrium states at finite temperature
- Zeeman-liouv to zeeman-hilb
- Sphten-liouv to zeeman-hilb
- Check the differences

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `equilibrium()`, `sphten2zeeman()`.
