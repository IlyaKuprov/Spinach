# examples/relaxation_theory/sle_nmr_dd_csa.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/sle_nmr_dd_csa.m`
- Signature: `sle_nmr_dd_csa()`
- Total lines: 82

## Purpose

15N-1H DD-CSA cross-correlation in a protein amide bond spin system using SLE formalism and Redfield relaxation theory. Calculation time: seconds

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: System specification; implemented by `sys.magnet=14.1`.
- Lines 25-26: Proximity cut-off; implemented by `sys.tols.prox_cutoff=4.0`.
- Lines 28-29: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 32-33: SLE housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 36-37: SLE parameters; implemented by `parameters.max_rank=10`.
- Lines 47-48: SLE simulation; implemented by `spectrum_sle=gridfree(spin_system,@slowpass,parameters,'nmr')`.
- Lines 50-51: SLE plotting; implemented by `kfigure(); subplot(1,2,1)`.
- Lines 55-56: BRW parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 61-62: BRW housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 65-66: BRW parameters; implemented by `parameters.spins={'15N'}`.
- Lines 73-74: BRW simulation; implemented by `spectrum_brw=liquid(spin_system,@slowpass,parameters,'nmr')`.
- Lines 76-77: BRW plotting; implemented by `subplot(1,2,2)`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 13: computes `sys.isotopes` using `sys.isotopes={'15N','1H'}`.
- Lines 14: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=[14.89 34.35 0.00`.
- Lines 17: computes `inter.zeeman.matrix{2}` using `inter.zeeman.matrix{2}=[30.75 1.31 0.00`.
- Lines 20: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(2)`.
- Lines 21: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=10`.
- Lines 22: computes `inter.coordinates` using `inter.coordinates={[-0.451455 -0.678015 0.000000]`.
- Lines 26: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 29: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 30: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 33: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 37: computes `parameters.max_rank` using `parameters.max_rank=10`.
- Lines 38: computes `parameters.tau_c` using `parameters.tau_c=5e-9`.
- Lines 39: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','15N')`.
- Lines 40: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','15N')`.
- Lines 41: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 42: computes `parameters.spins` using `parameters.spins={'15N'}`.
- Lines 43: computes `parameters.sweep` using `parameters.sweep=[-0.633e4 -0.629e4]`.

## Implementation structure

- 15N-1H DD-CSA cross-correlation in a protein amide bond
- spin system using SLE formalism and Redfield relaxation
- theory.
- Calculation time: seconds
- System specification
- Proximity cut-off
- Basis set
- SLE housekeeping
- SLE parameters
- SLE simulation
- SLE plotting
- BRW parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `gridfree()`, `kfigure()`, `subplot()`, `plot_1d()`, `ktitle()`, `kylabel()`, `liquid()`.
