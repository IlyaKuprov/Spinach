# tests/kernel/test_dynamic_voitlander.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_voitlander.m`
- Signature: `result=test_dynamic_voitlander()`
- Total lines: 191

## Purpose

Tests voitlander() on an isotropic one-electron field-swept line. Syntax: result=test_dynamic_voitlander()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Voitlander spherical-triangle integration\n')`.
- Lines 19-22: State the Voitlander target of the test; implemented by `result=new_test_result('kernel/dynamic_voitlander', 'Adaptive Voitlander triangle integral', 'voitlander() must integrate a constant isotropic transition as triangle are…`.
- Lines 24-25: Build an isotropic one-electron Hilbert-space system; implemented by `sys.magnet=1`.
- Lines 32-33: Set field-swept EPR parameters with a deliberately loose recursion gate; implemented by `parameters.spins={'E'}`.
- Lines 45-46: Get Zeeman, coupling, and microwave Hamiltonians; implemented by `[Ic,Qc]=hamiltonian(assume(spin_system,'labframe','couplings'))`.
- Lines 52-53: Find the isotropic transition at one orientation; implemented by `Hz=Iz+orientation(Qz,parameters.orientation)`.
- Lines 63-64: Define the positive octant spherical triangle; implemented by `r1=[1;0;0]`.
- Lines 68-69: Compute the first-level subdivision areas; implemented by `[r12,r23,r31]=sphtrsubd(r1,r2,r3)`.
- Lines 75-76: Keep all field points inside the internal six-width support window; implemented by `parameters.b_axis=tf+linspace(-2*tw,2*tw,parameters.npoints)`.
- Lines 78-82: Package the triangle vertices; implemented by `triangle=struct('xyz',{r1,r2,r3},'tf',{tf,tf,tf}, 'tm',{tm,tm,tm},'tw',{tw,tw,tw}, 'pd',{pd,pd,pd},'ti',{ti,ti,ti}, 'tj',{tj,tj,tj})`.
- Lines 84-85: Integrate the triangle using the production routine; implemented by `spec=voitlander(spin_system,parameters,triangle,Ic,Iz,Qc,Qz,Hmw)`.
- Lines 87-88: Build the analytic Lorentzian reference for a constant transition; implemented by `line_width=tw/2`.
- Lines 92-94: Check the returned spectrum against the constant-transition identity; implemented by `result=test_close(result,'voitlander constant line',spec,reference,1e-10,1e-12, 'subdivision of a constant transition must preserve total spherical-triangle area')`.
- Lines 99-100: Check that vertex-wise moment-Jacobian products are averaged directly; implemented by `tri_corr=triangle`.
- Lines 115-116: Check that singular field-sweep Jacobian vertices are rejected safely; implemented by `tri_sing=triangle`.
- Lines 131-132: Build a minimal two-level avoided crossing with two resonance roots; implemented by `spin_system.sys.output='hush'`.
- Lines 159-160: Get the two resonance roots and their generated branch labels; implemented by `tran=eigenfields(spin_system,parameters,Iz,Ic,Hmw)`.
- Lines 166-167: Compare stable and locally renumbered branch identities; implemented by `ti_ref=ti(high_root,:)`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/dynamic_voitlander', 'Adaptive Voitlander triangle integral', 'voitlander() must integrate a constant isotropic transition as triangle are…`.
- Lines 25: computes `sys.magnet` using `sys.magnet=1`.
- Lines 26: computes `sys.isotopes` using `sys.isotopes={'E'}`.
- Lines 27: computes `inter.zeeman.matrix` using `inter.zeeman.matrix={2.0023*eye(3)}`.
- Lines 28: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 29: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 30: computes `spin_system` using `spin_system=test_spin_system(sys,inter,bas)`.
- Lines 33: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 34: computes `parameters.mw_freq` using `parameters.mw_freq=9.5e9`.
- Lines 35: computes `parameters.fwhm` using `parameters.fwhm=2e-3`.
- Lines 36: computes `parameters.window` using `parameters.window=[0.33 0.35]`.
- Lines 37: computes `parameters.npoints` using `parameters.npoints=9`.
- Lines 38: computes `parameters.tm_tol` using `parameters.tm_tol=0`.
- Lines 39: computes `parameters.rspt_order` using `parameters.rspt_order=Inf`.
- Lines 40: computes `parameters.int_tol` using `parameters.int_tol=1e9`.
- Lines 41: computes `parameters.pp_tol` using `parameters.pp_tol=(parameters.window(2)-parameters.window(1))/(2*(parameters.npoints-1))`.
- Lines 42: computes `parameters.orientation` using `parameters.orientation=[0 0 0]`.
- Lines 43: computes `parameters.rho0` using `parameters.rho0=-state(spin_system,'Lz','E')`.

## Outputs

- result -regression test result with explanatory messages
- The test uses an isotropic spin-half electron, for which all triangle
- vertices and subdivision midpoints have the same transition field.

## Implementation structure

- Tests voitlander() on an isotropic one-electron field-swept line. Syntax:
- result=test_dynamic_voitlander()
- result -regression test result with explanatory messages
- The test uses an isotropic spin-half electron, for which all triangle
- vertices and subdivision midpoints have the same transition field.
- Announce the test target
- State the Voitlander target of the test
- Build an isotropic one-electron Hilbert-space system
- Set field-swept EPR parameters with a deliberately loose recursion gate
- Get Zeeman, coupling, and microwave Hamiltonians
- Find the isotropic transition at one orientation
- Define the positive octant spherical triangle

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `voitlander()`, `test_spin_system()`, `state()`, `hamiltonian()`, `assume()`, `orientation()`, `eigenfields()`, `test_true()`, `isscalar()`, `test_close()`, `sphtrsubd()`, `sphtarea()`, `all()`, `spec()`, `tri_corr()`.
