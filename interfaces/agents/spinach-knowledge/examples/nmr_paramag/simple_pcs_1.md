# examples/nmr_paramag/simple_pcs_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_paramag/simple_pcs_1.m`
- Signature: `simple_pcs_1()`
- Total lines: 50

## Purpose

Pseudocontact shift and Curie relaxation on a proton due to the presence of a point magnetic susceptibility centre. Calculation time: seconds

## Physical / mathematical content

- Paramagnetic NMR examples. These files work with pseudocontact shifts, paramagnetic relaxation, susceptibility tensors, and inverse problems for metal-site localisation or distributed electron density reconstruction.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: System specification; implemented by `sys.magnet=14.1`.
- Lines 14-15: Diamagnetic shifts and coordinates; implemented by `inter.zeeman.scalar={2.0}`.
- Lines 18-19: Magnetic susceptibility; implemented by `inter.suscept.chi={[0.0883 -0.0904 0.0822`.
- Lines 24-25: Relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 30-31: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 34-35: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 38-39: Spinach relaxation rates; implemented by `R=relaxation(spin_system)`.
- Lines 45-46: Summary; implemented by `disp(['R1 rate, Spinach: ' num2str(R1Sp)])`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 12: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 15: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={2.0}`.
- Lines 16: computes `inter.coordinates` using `inter.coordinates={[0.0 0.0 0.0]}`.
- Lines 19: computes `inter.suscept.chi` using `inter.suscept.chi={[0.0883 -0.0904 0.0822`.
- Lines 22: computes `inter.suscept.xyz` using `inter.suscept.xyz={[10.0 2.5 3.9]}`.
- Lines 25: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 26: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 27: computes `inter.rlx_keep` using `inter.rlx_keep='labframe'`.
- Lines 28: computes `inter.tau_c` using `inter.tau_c={10e-12}`.
- Lines 31: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 32: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 35: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 39: computes `R` using `R=relaxation(spin_system)`.
- Lines 40: computes `Lz` using `Lz=state(spin_system,'Lz','1H')`.
- Lines 41: computes `Lp` using `Lp=state(spin_system,'L+','1H')`.
- Lines 42: computes `R1Sp` using `R1Sp=-(Lz'*R*Lz)/(Lz'*Lz)`.
- Lines 43: computes `R2Sp` using `R2Sp=-(Lp'*R*Lp)/(Lp'*Lp)`.

## Implementation structure

- Pseudocontact shift and Curie relaxation on a proton due to
- the presence of a point magnetic susceptibility centre.
- Calculation time: seconds
- System specification
- Diamagnetic shifts and coordinates
- Magnetic susceptibility
- Relaxation theory parameters
- Basis set
- Spinach housekeeping
- Spinach relaxation rates
- Summary

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `relaxation()`, `state()`, `num2str()`.
