# examples/dnp_mas/solid_effect_mas_powder.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_mas/solid_effect_mas_powder.m`
- Signature: `solid_effect_mas_powder()`
- Total lines: 65

## Purpose

A MAS DNP simulation performed as described in Fred Mentink- Vigier's paper (Spinach rotation conventions are different): Steady state DNP simulation for a powder. Calculation time: minutes

## Physical / mathematical content

- MAS DNP examples. These files model microwave-driven electron-nuclear polarisation transfer under magic-angle spinning, combining rotor-synchronised anisotropic interactions, relaxation, microwave irradiation, and powder/rotor averaging.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Magnet field; implemented by `sys.magnet=9.403`.
- Lines 17-18: Spin specification; implemented by `sys.isotopes={'E','1H'}`.
- Lines 20-21: Interactions; implemented by `inter.zeeman.eigs{1}=[2.00614 2.00194 2.00988]`.
- Lines 28-29: Relaxation parameters; implemented by `inter.relaxation={'weizmann'}`.
- Lines 40-41: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 44-45: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 48-49: Stack generation parameters; implemented by `parameters.spins={'E'}`.
- Lines 60-61: Run the MAS DNP simulation; implemented by `answer=masdnp(spin_system,parameters)`.

### Key state/data transformations

- Lines 15: computes `sys.magnet` using `sys.magnet=9.403`.
- Lines 18: computes `sys.isotopes` using `sys.isotopes={'E','1H'}`.
- Lines 21: computes `inter.zeeman.eigs{1}` using `inter.zeeman.eigs{1}=[2.00614 2.00194 2.00988]`.
- Lines 22: computes `inter.zeeman.euler{1}` using `inter.zeeman.euler{1}=pi*[253.6 105.1 123.8]/180`.
- Lines 23: computes `inter.zeeman.eigs{2}` using `inter.zeeman.eigs{2}=[0.00 0.00 0.00]`.
- Lines 24: computes `inter.zeeman.euler{2}` using `inter.zeeman.euler{2}=[0.00 0.00 0.00]`.
- Lines 25-26: computes `inter.coordinates` using `inter.coordinates={[0.00 0.00 0.00], [0.00 0.00 3.00]}`.
- Lines 29: computes `inter.relaxation` using `inter.relaxation={'weizmann'}`.
- Lines 30: computes `inter.weiz_r1e` using `inter.weiz_r1e=1/0.3e-3`.
- Lines 31: computes `inter.weiz_r1n` using `inter.weiz_r1n=1/4.0`.
- Lines 32: computes `inter.weiz_r2e` using `inter.weiz_r2e=1/1.0e-6`.
- Lines 33: computes `inter.weiz_r2n` using `inter.weiz_r2n=1/0.2e-3`.
- Lines 34: computes `inter.weiz_r1d` using `inter.weiz_r1d=zeros(2,2)`.
- Lines 35: computes `inter.weiz_r2d` using `inter.weiz_r2d=zeros(2,2)`.
- Lines 36: computes `inter.temperature` using `inter.temperature=100`.
- Lines 37: computes `inter.equilibrium` using `inter.equilibrium='dibari'`.
- Lines 38: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 41: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.

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
- Stack generation parameters
- Run the MAS DNP simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `masdnp()`, `num2str()`.
