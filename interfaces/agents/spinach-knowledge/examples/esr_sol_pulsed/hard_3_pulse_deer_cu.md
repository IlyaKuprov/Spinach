# examples/esr_sol_pulsed/hard_3_pulse_deer_cu.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/hard_3_pulse_deer_cu.m`
- Signature: `hard_3_pulse_deer_cu()`
- Total lines: 72

## Purpose

Three-pulse DEER on a Cu(II)-NO two electron system at X-band. The numerical calculation is done by brute-force time propaga- tion and numerical powder averaging in Liouville space, inclu- ding g-factor orientation effects on the dipolar coupling. The analytical calculation is done for isotropic parts of the electron g-factors. Calculation time: seconds

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Spin system parameters; implemented by `sys.magnet=0.33`.
- Lines 25-26: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 29-30: Disable trajectory level SSR algorithms; implemented by `sys.disable={'trajlevel'}`.
- Lines 32-33: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 36-37: Sequence parameters; implemented by `parameters.rho0=state(spin_system,'Lz','E')`.
- Lines 47-48: Build the time axis; implemented by `time_axis=linspace(0,parameters.stepsize*parameters.nsteps,parameters.nsteps+1)`.
- Lines 50-51: Simulation (numerical); implemented by `deer_num=powder(spin_system,@deer_3p_hard_deer,parameters,'deer')`.
- Lines 53-54: Simulation (analytical); implemented by `D=xyz2dd(inter.coordinates{1},inter.coordinates{2},sys.isotopes{1},sys.isotopes{2})`.
- Lines 58-59: Plotting (numerical); implemented by `kfigure(); subplot(1,2,1)`.
- Lines 65-66: Plotting (analytical); implemented by `subplot(1,2,2); plot(1e6*time_axis,deer_anl)`.

### Key state/data transformations

- Lines 18: computes `sys.magnet` using `sys.magnet=0.33`.
- Lines 19: computes `sys.isotopes` using `sys.isotopes={'E','E'}`.
- Lines 20: computes `inter.zeeman.eigs` using `inter.zeeman.eigs={[2.056, 2.056, 2.205]`.
- Lines 22: computes `inter.zeeman.euler` using `inter.zeeman.euler={[0 0 0]; [0 0 0]}`.
- Lines 23: computes `inter.coordinates` using `inter.coordinates={[0 0 0]; [20 0 0]}`.
- Lines 26: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 27: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 30: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 33: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 37: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lz','E')`.
- Lines 38: computes `parameters.coil_prob` using `parameters.coil_prob=state(spin_system,{'L-'},{1})`.
- Lines 39: computes `parameters.stepsize` using `parameters.stepsize=1e-8`.
- Lines 40: computes `parameters.nsteps` using `parameters.nsteps=50`.
- Lines 41: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 42: computes `parameters.ex_prob` using `parameters.ex_prob=operator(spin_system,{'Lx'},{1})`.
- Lines 43: computes `parameters.ex_pump` using `parameters.ex_pump=operator(spin_system,{'Lx'},{2})`.
- Lines 44: computes `parameters.output` using `parameters.output='brief'`.
- Lines 45: computes `parameters.grid` using `parameters.grid='rep_2ang_1600pts_sph'`.

## Implementation structure

- Three-pulse DEER on a Cu(II)-NO two electron system at X-band.
- The numerical calculation is done by brute-force time propaga-
- tion and numerical powder averaging in Liouville space, inclu-
- ding g-factor orientation effects on the dipolar coupling.
- The analytical calculation is done for isotropic parts of the
- electron g-factors.
- Calculation time: seconds
- Spin system parameters
- Basis set
- Disable trajectory level SSR algorithms
- Spinach housekeeping
- Sequence parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `powder()`, `xyz2dd()`, `deer_analyt()`, `kfigure()`, `subplot()`, `kxlabel()`, `ktitle()`.
