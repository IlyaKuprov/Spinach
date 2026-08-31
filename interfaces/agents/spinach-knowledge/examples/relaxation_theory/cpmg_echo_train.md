# examples/relaxation_theory/cpmg_echo_train.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/cpmg_echo_train.m`
- Signature: `cpmg_echo_train()`
- Total lines: 55

## Purpose

CPMG echo train in a powder. Calculation time: seconds

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: System specification; implemented by `sys.magnet=14.1`.
- Lines 15-16: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 19-20: Relaxation theory; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 26-27: Algorithmic options; implemented by `sys.disable={'trajlevel'}`.
- Lines 29-30: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 33-34: Experiment setup; implemented by `parameters.grid='rep_2ang_200pts_sph'`.
- Lines 44-45: Simulation; implemented by `fid=powder(spin_system,@cpmg,parameters,'nmr')`.
- Lines 47-48: Plotting; implemented by `kfigure(); scale_figure([1.00 0.65])`.

### Key state/data transformations

- Lines 10: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 11: computes `sys.isotopes` using `sys.isotopes={'1H','1H'}`.
- Lines 12: computes `inter.zeeman.eigs` using `inter.zeeman.eigs={[-2 -2 4]-5,[-1 -3 4]+5}`.
- Lines 13: computes `inter.zeeman.euler` using `inter.zeeman.euler={[0 0 0],[0 0 0]}`.
- Lines 16: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 17: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 20: computes `inter.relaxation` using `inter.relaxation={'t1_t2'}`.
- Lines 21: computes `inter.r1_rates` using `inter.r1_rates={ 50.0 50.0}`.
- Lines 22: computes `inter.r2_rates` using `inter.r2_rates={150.0 150.0}`.
- Lines 23: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 24: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 27: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 30: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 34: computes `parameters.grid` using `parameters.grid='rep_2ang_200pts_sph'`.
- Lines 35: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 36: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','1H')`.
- Lines 37: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','1H')`.
- Lines 38: computes `parameters.pulse_op` using `parameters.pulse_op=operator(spin_system,'Lx','1H')`.

## Implementation structure

- CPMG echo train in a powder.
- Calculation time: seconds
- System specification
- Basis set
- Relaxation theory
- Algorithmic options
- Spinach housekeeping
- Experiment setup
- Simulation
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `powder()`, `kfigure()`, `scale_figure()`, `kxlabel()`, `kylabel()`.
