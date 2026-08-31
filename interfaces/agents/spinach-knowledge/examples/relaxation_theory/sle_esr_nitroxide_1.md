# examples/relaxation_theory/sle_esr_nitroxide_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/sle_esr_nitroxide_1.m`
- Signature: `sle_esr_nitroxide_1()`
- Total lines: 79

## Purpose

Comparison between nitroxide simulation using SLE formalism and Redfield relaxation theory. Calculation time: seconds

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Spin system properties; implemented by `options.no_xyz=1`.
- Lines 14-15: Magnet induction; implemented by `sys.magnet=3.5`.
- Lines 17-18: Proximity cut-off; implemented by `sys.tols.prox_cutoff=4.0`.
- Lines 20-21: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 24-25: SLE housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 28-29: SLE parameters; implemented by `parameters.max_rank=10`.
- Lines 41-42: SLE simulation; implemented by `spectrum_sle=gridfree(spin_system,@slowpass,parameters,'esr')`.
- Lines 45-46: SLE plotting; implemented by `kfigure(); subplot(1,2,1)`.
- Lines 50-51: BRW parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 56-57: BRW housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 60-61: BRW parameters; implemented by `parameters.spins={'E'}`.
- Lines 69-70: BRW simulation; implemented by `spectrum_brw=liquid(spin_system,@slowpass,parameters,'esr')`.
- Lines 73-74: BRW plotting; implemented by `subplot(1,2,2)`.

### Key state/data transformations

- Lines 11: computes `options.no_xyz` using `options.no_xyz=1`.
- Lines 12-13: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/nitroxide.log'), {{'E','E'},{'N','14N'}},[0 0],options)`.
- Lines 15: computes `sys.magnet` using `sys.magnet=3.5`.
- Lines 18: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 21: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 22: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 25: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 29: computes `parameters.max_rank` using `parameters.max_rank=10`.
- Lines 30: computes `parameters.tau_c` using `parameters.tau_c=5e-11`.
- Lines 31: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','E')`.
- Lines 32: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','E')`.
- Lines 33: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 34: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 35: computes `parameters.sweep` using `parameters.sweep=[-3e8 -1e8]`.
- Lines 36: computes `parameters.npoints` using `parameters.npoints=1024`.
- Lines 37: computes `parameters.zerofill` using `parameters.zerofill=1024`.
- Lines 38: computes `parameters.axis_units` using `parameters.axis_units='GHz-labframe'`.
- Lines 39: computes `parameters.invert_axis` using `parameters.invert_axis=1`.

## Implementation structure

- Comparison between nitroxide simulation using SLE formalism
- and Redfield relaxation theory.
- Calculation time: seconds
- Spin system properties
- Magnet induction
- Proximity cut-off
- Basis set
- SLE housekeeping
- SLE parameters
- SLE simulation
- SLE plotting
- BRW parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `state()`, `gridfree()`, `fdvec()`, `kfigure()`, `subplot()`, `plot_1d()`, `ktitle()`, `liquid()`.
