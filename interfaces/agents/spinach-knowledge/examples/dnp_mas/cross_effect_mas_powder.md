# examples/dnp_mas/cross_effect_mas_powder.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_mas/cross_effect_mas_powder.m`
- Signature: `cross_effect_mas_powder()`
- Total lines: 69

## Purpose

A MAS DNP simulation performed as described in Fred Mentink- Vigier's paper (Spinach rotation conventions are different): Steady state DNP simulation for a powder. Calculation time: minutes

## Physical / mathematical content

- MAS DNP examples. These files model microwave-driven electron-nuclear polarisation transfer under magic-angle spinning, combining rotor-synchronised anisotropic interactions, relaxation, microwave irradiation, and powder/rotor averaging.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Magnet field; implemented by `sys.magnet=9.394`.
- Lines 17-18: Spin specification; implemented by `sys.isotopes={'E','E','1H'}`.
- Lines 20-21: Interactions; implemented by `inter.zeeman.eigs{1}=[2.0094 2.0060 2.0017]`.
- Lines 34-35: Relaxation parameters; implemented by `inter.relaxation={'nottingham'}`.
- Lines 44-45: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 48-49: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 52-53: Experiment parameters; implemented by `parameters.spins={'E'}`.
- Lines 64-65: Run the MAS DNP simulation; implemented by `answer=masdnp(spin_system,parameters)`.

### Key state/data transformations

- Lines 15: computes `sys.magnet` using `sys.magnet=9.394`.
- Lines 18: computes `sys.isotopes` using `sys.isotopes={'E','E','1H'}`.
- Lines 21: computes `inter.zeeman.eigs{1}` using `inter.zeeman.eigs{1}=[2.0094 2.0060 2.0017]`.
- Lines 22: computes `inter.zeeman.euler{1}` using `inter.zeeman.euler{1}=[0.00 0.00 0.00]`.
- Lines 23: computes `inter.zeeman.eigs{2}` using `inter.zeeman.eigs{2}=[2.0094 2.0060 2.0017]`.
- Lines 24: computes `inter.zeeman.euler{2}` using `inter.zeeman.euler{2}=pi*[107 108 124]/180`.
- Lines 25: computes `inter.zeeman.eigs{3}` using `inter.zeeman.eigs{3}=[0.00 0.00 0.00]`.
- Lines 26: computes `inter.zeeman.euler{3}` using `inter.zeeman.euler{3}=[0.00 0.00 0.00]`.
- Lines 27: computes `inter.coupling.eigs` using `inter.coupling.eigs=cell(3,3)`.
- Lines 28: computes `inter.coupling.euler` using `inter.coupling.euler=cell(3,3)`.
- Lines 29: computes `inter.coupling.eigs{1,2}` using `inter.coupling.eigs{1,2}=[23.0e6 -11.5e6 -11.5e6]`.
- Lines 30: computes `inter.coupling.euler{1,2}` using `inter.coupling.euler{1,2}=pi*[0.00 135 0.00]/180`.
- Lines 31: computes `inter.coupling.eigs{1,3}` using `inter.coupling.eigs{1,3}=[1.5e6 -0.75e6 -0.75e6]`.
- Lines 32: computes `inter.coupling.euler{1,3}` using `inter.coupling.euler{1,3}=[0.00 0.00 0.00]`.
- Lines 35: computes `inter.relaxation` using `inter.relaxation={'nottingham'}`.
- Lines 36: computes `inter.nott_r1e` using `inter.nott_r1e=1/0.3e-3`.
- Lines 37: computes `inter.nott_r1n` using `inter.nott_r1n=1/4.0`.
- Lines 38: computes `inter.nott_r2e` using `inter.nott_r2e=1/1.0e-6`.

## Implementation structure

- A MAS DNP simulation performed as described in Fred Mentink-
- Vigier's paper (Spinach rotation conventions are different):
- Steady state DNP simulation for a powder.
- Calculation time: minutes
- Magnet field
- Spin specification
- Interactions
- Relaxation parameters
- Basis set
- Spinach housekeeping
- Experiment parameters
- Run the MAS DNP simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `masdnp()`, `num2str()`.
