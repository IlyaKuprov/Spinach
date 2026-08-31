# examples/fundamentals/state_tests/thermal_equilibrium_4.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/state_tests/thermal_equilibrium_4.m`
- Signature: `thermal_equilibrium_4()`
- Total lines: 72

## Purpose

Test of the invariance of the thermal equilibrium state under the thermalised relaxation superoperator.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 8-9: Magnet field; implemented by `sys.magnet=9.4`.
- Lines 11-12: Isotopes; implemented by `sys.isotopes={'19F','19F','19F','19F'}`.
- Lines 14-15: Chemical shifts; implemented by `inter.zeeman.scalar={-120.5380 -133.9429 -129.3169 -129.5320}`.
- Lines 17-18: J-couplings; implemented by `inter.coupling.scalar=cell(4,4)`.
- Lines 26-27: Relaxation theory parameters; implemented by `inter.relaxation={'damp'}`.
- Lines 33-34: Formalisms to test; implemented by `formalisms={'sphten-liouv','zeeman-liouv'}`.
- Lines 36-37: Loop over formalisms; implemented by `for n=1:numel(formalisms)`.
- Lines 39-40: Basis set; implemented by `bas.formalism=formalisms{n}`.
- Lines 43-44: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 47-48: Isotropic thermal equilibrium; implemented by `rho_eq=equilibrium(spin_system)`.
- Lines 50-51: Get relaxation superoperator; implemented by `R=relaxation(spin_system)`.
- Lines 53-54: Thermalise using Dibari-Levitt method; implemented by `HLSPS=hamiltonian(assume(spin_system,'labframe'),'left')`.
- Lines 57-58: Check the norm of the action; implemented by `if norm(Rt*rho_eq,2)>1e-9, error('R*rho_eq action test FAILED.'); end`.
- Lines 60-61: Thermalise using inhomogeneous master equation; implemented by `Rt=thermalize(spin_system,R,[],[],rho_eq,'IME')`.
- Lines 68-69: Report good test result; implemented by `disp('R*rho_eq action test PASSED.')`.

### Control flow inferred from the code

- Line 37: `for` loop over `n=1:numel(formalisms)`.
- Line 58: conditional branch on `norm(Rt*rho_eq,2)>1e-9, error('R*rho_eq action test FAILED.'); end`.
- Line 64: conditional branch on `norm(Rt*rho_eq,2)>1e-9, error('R*rho_eq action test FAILED.'); end`.

### Key state/data transformations

- Lines 9: computes `sys.magnet` using `sys.magnet=9.4`.
- Lines 12: computes `sys.isotopes` using `sys.isotopes={'19F','19F','19F','19F'}`.
- Lines 15: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={-120.5380 -133.9429 -129.3169 -129.5320}`.
- Lines 18: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(4,4)`.
- Lines 19: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}= 271.2924`.
- Lines 20: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}= 271.2924`.
- Lines 21: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}= 0.5401`.
- Lines 22: computes `inter.coupling.scalar{1,4}` using `inter.coupling.scalar{1,4}= -25.9884`.
- Lines 23: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}= 9.9625`.
- Lines 24: computes `inter.coupling.scalar{2,4}` using `inter.coupling.scalar{2,4}= -40.7675`.
- Lines 27: computes `inter.relaxation` using `inter.relaxation={'damp'}`.
- Lines 28: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 29: computes `inter.temperature` using `inter.temperature=40`.
- Lines 30: computes `inter.damp_rate` using `inter.damp_rate=5.0`.
- Lines 31: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 34: computes `formalisms` using `formalisms={'sphten-liouv','zeeman-liouv'}`.
- Lines 40: computes `bas.formalism` using `bas.formalism=formalisms{n}`.
- Lines 41: computes `bas.approximation` using `bas.approximation='none'`.

## Implementation structure

- Test of the invariance of the thermal equilibrium state under the
- thermalised relaxation superoperator.
- Magnet field
- Isotopes
- Chemical shifts
- J-couplings
- Relaxation theory parameters
- Formalisms to test
- Loop over formalisms
- Basis set
- Spinach housekeeping
- Isotropic thermal equilibrium

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `equilibrium()`, `relaxation()`, `hamiltonian()`, `assume()`, `thermalize()`.
