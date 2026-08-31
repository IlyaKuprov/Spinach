# examples/esr_sol_pulsed/ridme_cu_nitroxide.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/ridme_cu_nitroxide.m`
- Signature: `ridme_cu_nitroxide()`
- Total lines: 90

## Purpose

RIDME on a Cu(II)-NO two electron system at Q-band. The numerical calculation is done by brute-force time propaga- tion and numerical powder averaging in Liouville space, inclu- ding g-factor orientation effects on the dipolar coupling. Relaxation is included as extended T1/T2 theory. The analytical calculation is done for isotropic parts of the electron g-factors. Calculation time: seconds.

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 19-20: Spin system parameters; implemented by `sys.magnet=1.249`.
- Lines 27-28: Relaxation theory; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 34-35: Formalism; implemented by `bas.formalism='sphten-liouv'`.
- Lines 38-39: Disable trajectory level SSR algorithms; implemented by `sys.disable={'trajlevel'}`.
- Lines 41-42: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 45-46: Sequence parameters; implemented by `parameters.rho0=state(spin_system,'Lz','E')`.
- Lines 54-55: Simulation; implemented by `answer=powder(spin_system,@ridme,parameters,'esr')`.
- Lines 57-59: Phase cycle processing; implemented by `answer.ridme_trace.real=answer.pxpxpx.real+answer.pypypx.real+ answer.mxmxpx.real+answer.mymypx.real`.
- Lines 63-66: Axis ticks and figure frame; implemented by `x_axis=linspace((-parameters.stepsize*parameters.nsteps(1)), (parameters.stepsize*parameters.nsteps(2)), parameters.nsteps(1)+parameters.nsteps(2)+1)`.
- Lines 69-70: Real parts; implemented by `subplot(1,2,1); hold on; kgrid; box on`.
- Lines 79-80: Imaginary parts; implemented by `subplot(1,2,2); hold on; kgrid; box on`.

### Key state/data transformations

- Lines 20: computes `sys.magnet` using `sys.magnet=1.249`.
- Lines 21: computes `sys.isotopes` using `sys.isotopes={'E','E'}`.
- Lines 22: computes `inter.zeeman.eigs` using `inter.zeeman.eigs={[2.056, 2.056, 2.205]`.
- Lines 24: computes `inter.zeeman.euler` using `inter.zeeman.euler={[0 0 0]; [0 0 0]}`.
- Lines 25: computes `inter.coordinates` using `inter.coordinates={[43 0 0]; [0 0 0]}`.
- Lines 28: computes `inter.relaxation` using `inter.relaxation={'t1_t2'}`.
- Lines 29: computes `inter.r1_rates` using `inter.r1_rates={1/(35e-6) 1/(2e-3)}`.
- Lines 30: computes `inter.r2_rates` using `inter.r2_rates={1/(1.5e-6) 1/(1.3e-6)}`.
- Lines 31: computes `inter.rlx_keep` using `inter.rlx_keep='labframe'`.
- Lines 32: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 35: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 36: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 39: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 42: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 46: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lz','E')`.
- Lines 47: computes `parameters.probe_spin` using `parameters.probe_spin=2`.
- Lines 48: computes `parameters.stepsize` using `parameters.stepsize=16e-9`.
- Lines 49: computes `parameters.nsteps` using `parameters.nsteps=[25 188]`.

## Implementation structure

- RIDME on a Cu(II)-NO two electron system at Q-band.
- The numerical calculation is done by brute-force time propaga-
- tion and numerical powder averaging in Liouville space, inclu-
- ding g-factor orientation effects on the dipolar coupling.
- Relaxation is included as extended T1/T2 theory.
- The analytical calculation is done for isotropic parts of the
- electron g-factors.
- Calculation time: seconds.
- Spin system parameters
- Relaxation theory
- Formalism
- Disable trajectory level SSR algorithms

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `powder()`, `kfigure()`, `scale_figure()`, `subplot()`, `kxlabel()`, `ktitle()`, `klegend()`.
