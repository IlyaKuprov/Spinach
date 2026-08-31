# examples/esr_liq_pulsed/rapidscan_nitroxide.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_liq_pulsed/rapidscan_nitroxide.m`
- Signature: `rapidscan_nitroxide()`
- Total lines: 55

## Purpose

Rapid scan ESR spectrum of a nitroxide radical. Calculation time: seconds

## Physical / mathematical content

- Liquid-state ESR examples. The dominant physics is electron Zeeman interaction, hyperfine coupling, relaxation broadening, and pulse-acquire or ENDOR-type detection in fast tumbling systems.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Centre field; implemented by `sys.magnet=3.5`.
- Lines 12-13: Spin system properties; implemented by `sys.isotopes={'14N','E'}`.
- Lines 27-28: Simulation parameters; implemented by `bas.formalism='sphten-liouv'`.
- Lines 36-37: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 40-41: Experiment parameters; implemented by `parameters.mw_pwr=2*pi*1e3`.
- Lines 46-47: Run the experiment; implemented by `[x_axis,spectrum]=rapidscan(spin_system,parameters)`.
- Lines 49-50: Plot the result; implemented by `kfigure(); plot(x_axis,real(spectrum)); kgrid`.

### Key state/data transformations

- Lines 10: computes `sys.magnet` using `sys.magnet=3.5`.
- Lines 13: computes `sys.isotopes` using `sys.isotopes={'14N','E'}`.
- Lines 14: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=zeros(3)`.
- Lines 15: computes `inter.zeeman.matrix{2}` using `inter.zeeman.matrix{2}=[2.0104 0.0000 0.0001`.
- Lines 18: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=[]`.
- Lines 19: computes `inter.coupling.matrix{2,2}` using `inter.coupling.matrix{2,2}=[]`.
- Lines 20: computes `inter.coupling.matrix{1,2}` using `inter.coupling.matrix{1,2}=[0.6178 0 0.3161`.
- Lines 23: computes `inter.coupling.matrix{2,1}` using `inter.coupling.matrix{2,1}=[0.6178 0 0.3161`.
- Lines 28: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 29: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 30: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 31: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 32: computes `inter.equilibrium` using `inter.equilibrium='dibari'`.
- Lines 33: computes `inter.temperature` using `inter.temperature=100`.
- Lines 34: computes `inter.tau_c` using `inter.tau_c={2e-11}`.
- Lines 37: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 41: computes `parameters.mw_pwr` using `parameters.mw_pwr=2*pi*1e3`.
- Lines 42: computes `parameters.sweep` using `parameters.sweep=[-0.011 -0.003]`.

## Implementation structure

- Rapid scan ESR spectrum of a nitroxide radical.
- Calculation time: seconds
- Centre field
- Spin system properties
- Simulation parameters
- Spinach housekeeping
- Experiment parameters
- Run the experiment
- Plot the result

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `rapidscan()`, `kfigure()`, `kxlabel()`, `kylabel()`.
