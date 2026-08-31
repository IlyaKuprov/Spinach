# examples/fundamentals/derivative_tests/dirdiff_test_system.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/derivative_tests/dirdiff_test_system.m`
- Signature: `[spin_system,Sx,Sy,Sz,Lx,Ly,H]=dirdiff_test_system(formalism)`
- Total lines: 80

## Purpose

Spin system generator for directional derivative tests. Syntax: [spin_system,Sx,Sy,Sz,Lx,Ly,H]=dirdiff_test_system(formalism)

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Check consistency; implemented by `grumble(formalism)`.
- Lines 12-13: Select system size; implemented by `switch formalism`.
- Lines 17-18: Keep the original large Liouville-space test; implemented by `n_spins=100`.
- Lines 22-23: Use a compact system for full Zeeman formalisms; implemented by `n_spins=2`.
- Lines 27-28: Set the magnetic field; implemented by `sys.magnet=28.18`.
- Lines 31-33: Put non-interacting spins at equal intervals within the [-100,+100] ppm chemical shift range; implemented by `sys.isotopes=cell(n_spins,1)`.
- Lines 39-40: Select the requested basis set; implemented by `bas.formalism=formalism`.
- Lines 45-46: Keep complete single-spin terms only; implemented by `bas.approximation='IK-2'`.
- Lines 52-53: Keep the full Zeeman basis; implemented by `bas.approximation='none'`.
- Lines 57-58: Run Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 61-62: Set up spin states; implemented by `Sx=state(spin_system,'Lx','13C'); Sx=Sx/norm(full(Sx),2)`.
- Lines 66-67: Get control operators; implemented by `Lx=operator(spin_system,'Lx','13C')`.
- Lines 70-71: Get the drift Hamiltonian; implemented by `H=hamiltonian(assume(spin_system,'nmr'))`.

### Control flow inferred from the code

- Line 13: dispatches on `formalism`; cases `'sphten-liouv'`, `{'zeeman-liouv','zeeman-hilb'}`, `'sphten-liouv'`, `{'zeeman-liouv','zeeman-hilb'}`.
- Line 34: `for` loop over `n=1:n_spins`.
- Line 41: dispatches on `formalism`; cases `'sphten-liouv'`, `{'zeeman-liouv','zeeman-hilb'}`.

### Key state/data transformations

- Lines 18: computes `n_spins` using `n_spins=100`.
- Lines 28: computes `sys.magnet` using `sys.magnet=28.18`.
- Lines 29: computes `sys.output` using `sys.output='hush'`.
- Lines 33: computes `sys.isotopes` using `sys.isotopes=cell(n_spins,1)`.
- Lines 35: computes `sys.isotopes{n}` using `sys.isotopes{n}='13C'`.
- Lines 37: computes `inter.zeeman.scalar` using `inter.zeeman.scalar=num2cell(linspace(-100,100,n_spins))`.
- Lines 40: computes `bas.formalism` using `bas.formalism=formalism`.
- Lines 46: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 47: computes `bas.space_level` using `bas.space_level=1`.
- Lines 48: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 58: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 62: computes `Sx` using `Sx=state(spin_system,'Lx','13C'); Sx=Sx/norm(full(Sx),2)`.
- Lines 63: computes `Sy` using `Sy=state(spin_system,'Ly','13C'); Sy=Sy/norm(full(Sy),2)`.
- Lines 64: computes `Sz` using `Sz=state(spin_system,'Lz','13C'); Sz=Sz/norm(full(Sz),2)`.
- Lines 67: computes `Lx` using `Lx=operator(spin_system,'Lx','13C')`.
- Lines 68: computes `Ly` using `Ly=operator(spin_system,'Ly','13C')`.
- Lines 71: computes `H` using `H=hamiltonian(assume(spin_system,'nmr'))`.

### Local helper functions

- Line 76: `grumble()` — `function grumble(formalism)`.
  - Representative operation: `if (~ischar(formalism))||(~ismember(formalism,{'sphten-liouv','zeeman-liouv','zeeman-hilb'}))`.
  - Representative operation: `error('formalism must be ''sphten-liouv'', ''zeeman-liouv'', or ''zeeman-hilb''.')`.

## Implementation structure

- Spin system generator for directional derivative tests. Syntax:
- [spin_system,Sx,Sy,Sz,Lx,Ly,H]=dirdiff_test_system(formalism)
- Check consistency
- Select system size
- Keep the original large Liouville-space test
- Use a compact system for full Zeeman formalisms
- Set the magnetic field
- Put non-interacting spins at equal intervals
- within the [-100,+100] ppm chemical shift range
- Select the requested basis set
- Keep complete single-spin terms only
- Keep the full Zeeman basis

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `num2cell()`, `create()`, `basis()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `ischar()`, `ismember()`.
