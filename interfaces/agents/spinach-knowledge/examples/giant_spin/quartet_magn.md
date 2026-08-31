# examples/giant_spin/quartet_magn.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/quartet_magn.m`
- Signature: `quartet_magn()`
- Total lines: 52

## Purpose

Sample magnetisation during a finite-speed magnetic field sweep for a spin-3/2 particle with a zero-field splitting. Calculation time: seconds

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
- Lines 28-29: Temperature; implemented by `inter.temperature=1.0`.
- Lines 31-32: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 35-36: Experiment parameters; implemented by `parameters.fields=[0 1]`.
- Lines 42-43: Run the field scan; implemented by `[fields,z_magn]=fieldscan_magn(spin_system,parameters)`.
- Lines 45-46: Plot the results; implemented by `kfigure(); plot(fields,z_magn); kgrid`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=1.0`.
- Lines 15: computes `sys.isotopes` using `sys.isotopes={'E4'}`.
- Lines 18: computes `inter.zeeman.matrix` using `inter.zeeman.matrix={[2 0 0; 0 2 0; 0 0 2]}`.
- Lines 21: computes `D` using `D=icm2hz(-0.5); E=0.3*D`.
- Lines 22: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=zfs2mat(D,E,0,0,0)`.
- Lines 25: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 26: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 29: computes `inter.temperature` using `inter.temperature=1.0`.
- Lines 32: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 36: computes `parameters.fields` using `parameters.fields=[0 1]`.
- Lines 37: computes `parameters.npoints` using `parameters.npoints=1000`.
- Lines 38: computes `parameters.sweep_time` using `parameters.sweep_time=1e-9`.
- Lines 39: computes `parameters.orientation` using `parameters.orientation=[0 0 0]`.
- Lines 40: computes `parameters.nstates` using `parameters.nstates=4`.
- Lines 43: computes `[fields,z_magn]` using `[fields,z_magn]=fieldscan_magn(spin_system,parameters)`.

## Implementation structure

- Sample magnetisation during a finite-speed magnetic field
- sweep for a spin-3/2 particle with a zero-field splitting.
- Calculation time: seconds
- This must be set to 1 Tesla
- Particle
- Zeeman tensor
- Zero-field splitting
- Formalism and basis set
- Temperature
- Spinach housekeeping
- Experiment parameters
- Run the field scan

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `icm2hz()`, `zfs2mat()`, `create()`, `basis()`, `fieldscan_magn()`, `kfigure()`, `kxlabel()`, `kylabel()`.
