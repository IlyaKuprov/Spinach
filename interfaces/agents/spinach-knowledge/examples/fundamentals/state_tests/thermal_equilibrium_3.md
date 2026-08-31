# examples/fundamentals/state_tests/thermal_equilibrium_3.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/state_tests/thermal_equilibrium_3.m`
- Signature: `thermal_equilibrium_3()`
- Total lines: 73

## Purpose

Test of the thermal equilibrium functionality against the textbook expressions for the Boltzmann populations.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 8-9: X-band magnet; implemented by `sys.magnet=0.34`.
- Lines 11-12: Electron and two protons; implemented by `sys.isotopes={'E','1H','1H'}`.
- Lines 14-15: Zeeman interactions (g-tensor for trityl, ppm guess for 1H); implemented by `inter.zeeman.eigs={[2.00319 2.00319 2.00258],[0 0 5],[0 5 0]}`.
- Lines 18-19: Cartesian coordinates; implemented by `inter.coordinates={[0.000 0.000 0.000]`.
- Lines 23-24: Spin temperature; implemented by `inter.temperature=80`.
- Lines 26-27: Formalisms to test; implemented by `formalisms={'zeeman-hilb','zeeman-liouv','sphten-liouv'}`.
- Lines 29-30: Loop over formalisms; implemented by `for n=1:numel(formalisms)`.
- Lines 32-33: Formalism and basis set; implemented by `bas.formalism=formalisms{n}`.
- Lines 36-37: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 40-41: Isotropic thermal equilibrium; implemented by `rho_eq=equilibrium(spin_system)`.
- Lines 43-44: Detection states; implemented by `coil_a=state(spin_system,{'Lz'},{1})`.
- Lines 48-49: Get the expectation values from Spinach; implemented by `expt_a_spinach=trace(coil_a'*rho_eq)`.
- Lines 53-54: Get the expectation values from textbook; implemented by `[~,PE]=levelpop('E',sys.magnet,inter.temperature)`.
- Lines 60-63: Test for differences; implemented by `if (abs((expt_a_spinach-expt_a_textbook)/expt_a_textbook)>1e-3)|| (abs((expt_b_spinach-expt_b_textbook)/expt_b_textbook)>1e-3)|| (abs((expt_c_spinach-expt_c_textbook)/ex…`.
- Lines 69-70: Report success; implemented by `disp('Thermodynamic equilibrium tests PASSED.')`.

### Control flow inferred from the code

- Line 30: `for` loop over `n=1:numel(formalisms)`.
- Line 61: conditional branch on `(abs((expt_a_spinach-expt_a_textbook)/expt_a_textbook)>1e-3)||`.

### Key state/data transformations

- Lines 9: computes `sys.magnet` using `sys.magnet=0.34`.
- Lines 12: computes `sys.isotopes` using `sys.isotopes={'E','1H','1H'}`.
- Lines 15: computes `inter.zeeman.eigs` using `inter.zeeman.eigs={[2.00319 2.00319 2.00258],[0 0 5],[0 5 0]}`.
- Lines 16: computes `inter.zeeman.euler` using `inter.zeeman.euler=(pi/180)*{[0 10 0],[0 0 10],[100 0 0]}`.
- Lines 19: computes `inter.coordinates` using `inter.coordinates={[0.000 0.000 0.000]`.
- Lines 24: computes `inter.temperature` using `inter.temperature=80`.
- Lines 27: computes `formalisms` using `formalisms={'zeeman-hilb','zeeman-liouv','sphten-liouv'}`.
- Lines 33: computes `bas.formalism` using `bas.formalism=formalisms{n}`.
- Lines 34: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 37: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 41: computes `rho_eq` using `rho_eq=equilibrium(spin_system)`.
- Lines 44: computes `coil_a` using `coil_a=state(spin_system,{'Lz'},{1})`.
- Lines 45: computes `coil_b` using `coil_b=state(spin_system,{'Lz'},{2})`.
- Lines 46: computes `coil_c` using `coil_c=state(spin_system,{'Lz'},{3})`.
- Lines 49: computes `expt_a_spinach` using `expt_a_spinach=trace(coil_a'*rho_eq)`.
- Lines 50: computes `expt_b_spinach` using `expt_b_spinach=trace(coil_b'*rho_eq)`.
- Lines 51: computes `expt_c_spinach` using `expt_c_spinach=trace(coil_c'*rho_eq)`.
- Lines 54: computes `[~,PE]` using `[~,PE]=levelpop('E',sys.magnet,inter.temperature)`.

## Implementation structure

- Test of the thermal equilibrium functionality against the
- textbook expressions for the Boltzmann populations.
- X-band magnet
- Electron and two protons
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
- Cartesian coordinates
- Spin temperature
- Formalisms to test
- Loop over formalisms
- Formalism and basis set
- Spinach housekeeping
- Isotropic thermal equilibrium

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `equilibrium()`, `state()`, `levelpop()`.
