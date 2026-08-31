# examples/fundamentals/state_tests/thermal_equilibrium_5.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/state_tests/thermal_equilibrium_5.m`
- Signature: `thermal_equilibrium_5()`
- Total lines: 86

## Purpose

Cross-formalism test of state recovery towards the thermodynamic equilibrium.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 8-9: Magnet field; implemented by `sys.magnet=9.4`.
- Lines 11-12: Isotopes; implemented by `sys.isotopes={'19F','19F','19F','19F'}`.
- Lines 14-15: Chemical shifts; implemented by `inter.zeeman.scalar={-120.5380 -133.9429 -129.3169 -129.5320}`.
- Lines 17-18: J-couplings; implemented by `inter.coupling.scalar=cell(4,4)`.
- Lines 26-27: Relaxation theory parameters; implemented by `inter.relaxation={'damp'}`.
- Lines 32-33: Formalisms and methods to test; implemented by `formalisms={'sphten-liouv','zeeman-liouv'}`.
- Lines 36-37: Loop over formalisms; implemented by `for n=1:numel(formalisms)`.
- Lines 39-40: Loop over methods; implemented by `for k=1:numel(methods)`.
- Lines 42-43: Thermalisation method; implemented by `inter.equilibrium=methods{k}`.
- Lines 45-46: Basis set; implemented by `bas.formalism=formalisms{n}`.
- Lines 49-50: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 53-54: Get thermal equilibrium state; implemented by `rho_eq=equilibrium(spin_system)`.
- Lines 56-57: Flip the magnetisation over; implemented by `Fx=operator(spin_system,'Lx','19F')`.
- Lines 60-61: Get the evolution generator; implemented by `spin_system=assume(spin_system,'nmr')`.
- Lines 64-65: Get recovery trajectories; implemented by `coil=state(spin_system,'Lz','19F')`.
- Lines 72-73: Differences between formalisms and methods; implemented by `a=norm(traj{1,1}-traj{2,1},2)/norm(traj{1,1}+traj{2,1},2)`.
- Lines 78-79: Report the test result; implemented by `if max([a b c d])>1e-3`.

### Control flow inferred from the code

- Line 37: `for` loop over `n=1:numel(formalisms)`.
- Line 40: `for` loop over `k=1:numel(methods)`.
- Line 79: conditional branch on `max([a b c d])>1e-3`.

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
- Lines 28: computes `inter.temperature` using `inter.temperature=40`.
- Lines 29: computes `inter.damp_rate` using `inter.damp_rate=5.0`.
- Lines 30: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 33: computes `formalisms` using `formalisms={'sphten-liouv','zeeman-liouv'}`.
- Lines 34: computes `methods` using `methods={'dibari','IME'}; traj=cell(2,2)`.
- Lines 43: computes `inter.equilibrium` using `inter.equilibrium=methods{k}`.
- Lines 46: computes `bas.formalism` using `bas.formalism=formalisms{n}`.

## Implementation structure

- Cross-formalism test of state recovery towards the
- thermodynamic equilibrium.
- Magnet field
- Isotopes
- Chemical shifts
- J-couplings
- Relaxation theory parameters
- Formalisms and methods to test
- Loop over formalisms
- Loop over methods
- Thermalisation method
- Basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `equilibrium()`, `operator()`, `step()`, `assume()`, `hamiltonian()`, `relaxation()`, `state()`, `evolution()`.
