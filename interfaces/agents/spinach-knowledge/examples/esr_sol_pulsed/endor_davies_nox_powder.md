# examples/esr_sol_pulsed/endor_davies_nox_powder.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/endor_davies_nox_powder.m`
- Signature: `endor_davies_nox_powder()`
- Total lines: 79

## Purpose

Davies ENDOR simulation for a nitroxide radical. Soft pulses are simulated using Fokker-Planck formalism. This is a pain- fully slow brute-force time-domain simulation with explicit soft pulses and full account of the effect of the orientati- on selection using a large spherical averaging grid. Calculation time: hours.

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Isotopes; implemented by `sys.isotopes={'E','14N'}`.
- Lines 16-17: Magnet field; implemented by `sys.magnet=3.5`.
- Lines 19-20: Interactions; implemented by `inter.zeeman.matrix=cell(1,2)`.
- Lines 29-30: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 33-34: Relaxation theory; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 40-41: Disable trajectory-level SSR algorithms; implemented by `sys.disable={'trajlevel'}`.
- Lines 43-44: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 47-48: Sequence parameters; implemented by `parameters.spins={'E','14N'}`.
- Lines 56-57: Electron pulse parameters; implemented by `parameters.e_rnk=2`.
- Lines 63-64: Nucleus pulse parameters; implemented by `parameters.n_rnk=3`.
- Lines 70-71: Simulation; implemented by `answer=powder(spin_system,@endor_davies,parameters,'esr')`.
- Lines 73-74: Plotting; implemented by `kfigure(); plot(parameters.n_frq/1e6,real(answer))`.

### Key state/data transformations

- Lines 14: computes `sys.isotopes` using `sys.isotopes={'E','14N'}`.
- Lines 17: computes `sys.magnet` using `sys.magnet=3.5`.
- Lines 20: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1,2)`.
- Lines 21: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=[2.01045 0.00000 0.00000`.
- Lines 24: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(2,2)`.
- Lines 25: computes `inter.coupling.matrix{1,2}` using `inter.coupling.matrix{1,2}=[1.2356 0.0000 0.6322`.
- Lines 30: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 31: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 34: computes `inter.relaxation` using `inter.relaxation={'t1_t2'}`.
- Lines 35: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 36: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 37: computes `inter.r1_rates` using `inter.r1_rates={20e6 0.5e6}`.
- Lines 38: computes `inter.r2_rates` using `inter.r2_rates={20e6 0.5e6}`.
- Lines 41: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 44: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 48: computes `parameters.spins` using `parameters.spins={'E','14N'}`.
- Lines 49: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lz','E')`.
- Lines 50: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','E')`.

## Implementation structure

- Davies ENDOR simulation for a nitroxide radical. Soft pulses
- are simulated using Fokker-Planck formalism. This is a pain-
- fully slow brute-force time-domain simulation with explicit
- soft pulses and full account of the effect of the orientati-
- on selection using a large spherical averaging grid.
- Calculation time: hours.
- Isotopes
- Magnet field
- Interactions
- Basis set
- Relaxation theory
- Disable trajectory-level SSR algorithms

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `powder()`, `kfigure()`, `kylabel()`, `kxlabel()`.
