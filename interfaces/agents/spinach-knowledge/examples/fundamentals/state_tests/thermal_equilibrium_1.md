# examples/fundamentals/state_tests/thermal_equilibrium_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/state_tests/thermal_equilibrium_1.m`
- Signature: `thermal_equilibrium_1()`
- Total lines: 57

## Purpose

Observables at thermal equilibrium using the three formalisms supported by Spinach kernel, tested against known answers.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 8-9: Spin system parameters; implemented by `sys.magnet=14.1`.
- Lines 20-21: Preallocate the answers; implemented by `eq_mags_spinach=zeros(4,3,'like',1i)`.
- Lines 24-25: Get numerical equilibrium magnetisation; implemented by `formalisms={'zeeman-hilb','zeeman-liouv','sphten-liouv'}`.
- Lines 38-39: Get analytical equilibrium magnetisations; implemented by `[~,P]=levelpop('E8',sys.magnet,inter.temperature)`.
- Lines 49-50: Display the answers; implemented by `if any(abs(eq_mags_spinach-eq_mags_textbook)>1e-5,'all')`.

### Control flow inferred from the code

- Line 26: `for` loop over `n=1:numel(formalisms)`.
- Line 50: conditional branch on `any(abs(eq_mags_spinach-eq_mags_textbook)>1e-5,'all')`.

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
- Lines 21: computes `eq_mags_spinach` using `eq_mags_spinach=zeros(4,3,'like',1i)`.
- Lines 22: computes `eq_mags_textbook` using `eq_mags_textbook=zeros(4,1,'like',1i)`.
- Lines 25: computes `formalisms` using `formalisms={'zeeman-hilb','zeeman-liouv','sphten-liouv'}`.
- Lines 27: computes `bas.formalism` using `bas.formalism=formalisms{n}`.
- Lines 28: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 29: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 31: computes `rho` using `rho=equilibrium(spin_system)`.
- Lines 32: computes `eq_mags_spinach(1,n)` using `eq_mags_spinach(1,n)=trace(state(spin_system,'Lz','E8')'*rho)`.

## Implementation structure

- Observables at thermal equilibrium using the three formalisms
- supported by Spinach kernel, tested against known answers.
- Spin system parameters
- Preallocate the answers
- Get numerical equilibrium magnetisation
- Get analytical equilibrium magnetisations
- Display the answers

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `equilibrium()`, `eq_mags_spinach()`, `state()`, `levelpop()`, `eq_mags_textbook()`, `any()`.
