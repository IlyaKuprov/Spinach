# examples/giant_spin/quartet_levels.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/quartet_levels.m`
- Signature: `quartet_levels()`
- Total lines: 42

## Purpose

Energy levels magnetic field scan for a spin-3/2 particle with a zero-field splitting. Calculation time: seconds

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: This must be set to 1 Tesla; implemented by `sys.magnet=1.0`.
- Lines 14-15: Particle; implemented by `sys.isotopes={'E4'}`.
- Lines 17-18: Zeeman tensor; implemented by `inter.zeeman.matrix={[2 0 0; 0 2 0; 0 0 2]}`.
- Lines 20-21: Zero-field splitting; implemented by `D=icm2hz(-0.5); E=0.3*D`.
- Lines 24-25: Formalism and basis set; implemented by `bas.approximation='none'`.
- Lines 28-29: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 32-33: Experiment parameters; implemented by `parameters.fields=[0 1]`.
- Lines 38-39: Run the field scan; implemented by `fieldscan_enlev(spin_system,parameters)`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=1.0`.
- Lines 15: computes `sys.isotopes` using `sys.isotopes={'E4'}`.
- Lines 18: computes `inter.zeeman.matrix` using `inter.zeeman.matrix={[2 0 0; 0 2 0; 0 0 2]}`.
- Lines 21: computes `D` using `D=icm2hz(-0.5); E=0.3*D`.
- Lines 22: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=zfs2mat(D,E,0,0,0)`.
- Lines 25: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 26: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 29: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 33: computes `parameters.fields` using `parameters.fields=[0 1]`.
- Lines 34: computes `parameters.npoints` using `parameters.npoints=100`.
- Lines 35: computes `parameters.orientation` using `parameters.orientation=[0 0 0]`.
- Lines 36: computes `parameters.nstates` using `parameters.nstates=4`.

## Implementation structure

- Energy levels magnetic field scan for a spin-3/2 particle
- with a zero-field splitting.
- Calculation time: seconds
- This must be set to 1 Tesla
- Particle
- Zeeman tensor
- Zero-field splitting
- Formalism and basis set
- Spinach housekeeping
- Experiment parameters
- Run the field scan

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `icm2hz()`, `zfs2mat()`, `create()`, `basis()`, `fieldscan_enlev()`.
