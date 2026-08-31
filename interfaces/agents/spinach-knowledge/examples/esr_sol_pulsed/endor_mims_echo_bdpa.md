# examples/esr_sol_pulsed/endor_mims_echo_bdpa.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/endor_mims_echo_bdpa.m`
- Signature: `endor_mims_echo_bdpa()`
- Total lines: 57

## Purpose

Stimulated echo stage of the Mims ENDOR pulse sequence on BDPA. The nuclear pulse is not applied, this is echo dia- gnostics stage. The echo gets sharper when g-tensor aniso- tropy is increased. Run time: seconds.

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Isotopes; implemented by `sys.isotopes={'E','1H','1H'}`.
- Lines 16-17: Magnet field; implemented by `sys.magnet=3.35`.
- Lines 19-20: Interactions; implemented by `inter.zeeman.matrix=cell(1,3)`.
- Lines 32-33: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 36-37: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 40-41: Sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 48-49: Simulation; implemented by `answer=powder(spin_system,@endor_mims_echo,parameters,'esr')`.
- Lines 51-52: Plotting; implemented by `time_axis=linspace(0,2*parameters.tau,parameters.nsteps+1)`.

### Key state/data transformations

- Lines 14: computes `sys.isotopes` using `sys.isotopes={'E','1H','1H'}`.
- Lines 17: computes `sys.magnet` using `sys.magnet=3.35`.
- Lines 20: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1,3)`.
- Lines 21: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=[2.00263 0.00000 0.00000`.
- Lines 24: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(3,3)`.
- Lines 25: computes `inter.coupling.matrix{1,2}` using `inter.coupling.matrix{1,2}=[7.70 0.00 0.00`.
- Lines 28: computes `inter.coupling.matrix{1,3}` using `inter.coupling.matrix{1,3}=[1.00 0.00 0.00`.
- Lines 33: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 34: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 37: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 41: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 42: computes `parameters.electrons` using `parameters.electrons=1`.
- Lines 43: computes `parameters.tau` using `parameters.tau=200e-9`.
- Lines 44: computes `parameters.n_dur` using `parameters.n_dur=50e-6`.
- Lines 45: computes `parameters.grid` using `parameters.grid='rep_2ang_400pts_sph'`.
- Lines 46: computes `parameters.nsteps` using `parameters.nsteps=200`.
- Lines 49: computes `answer` using `answer=powder(spin_system,@endor_mims_echo,parameters,'esr')`.
- Lines 52: computes `time_axis` using `time_axis=linspace(0,2*parameters.tau,parameters.nsteps+1)`.

## Implementation structure

- Stimulated echo stage of the Mims ENDOR pulse sequence on
- BDPA. The nuclear pulse is not applied, this is echo dia-
- gnostics stage. The echo gets sharper when g-tensor aniso-
- tropy is increased.
- Run time: seconds.
- Isotopes
- Magnet field
- Interactions
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `powder()`, `kfigure()`, `kxlabel()`, `kylabel()`.
