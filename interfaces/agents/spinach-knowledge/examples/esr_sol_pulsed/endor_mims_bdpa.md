# examples/esr_sol_pulsed/endor_mims_bdpa.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/endor_mims_bdpa.m`
- Signature: `endor_mims_bdpa()`
- Total lines: 67

## Purpose

Mims ENDOR pulse sequence on BDPA with ideal electron pulses, reproducing Figure 10 from Calculation time: hours, much faster on GPU

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Isotopes; implemented by `sys.isotopes={'E','1H','1H'}`.
- Lines 16-17: Magnet field; implemented by `sys.magnet=3.35`.
- Lines 19-20: Interactions; implemented by `inter.zeeman.matrix=cell(1,3)`.
- Lines 32-33: Relaxation theory; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 39-40: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 43-44: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 47-48: Sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 58-59: Simulation; implemented by `answer=powder(spin_system,@endor_mims_ideal,parameters,'esr')`.
- Lines 61-62: Plotting; implemented by `kfigure(); plot(1e-6*parameters.n_frq,abs(answer))`.

### Key state/data transformations

- Lines 14: computes `sys.isotopes` using `sys.isotopes={'E','1H','1H'}`.
- Lines 17: computes `sys.magnet` using `sys.magnet=3.35`.
- Lines 20: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1,3)`.
- Lines 21: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=[2.00263 0.00000 0.00000`.
- Lines 24: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(3,3)`.
- Lines 25: computes `inter.coupling.matrix{1,2}` using `inter.coupling.matrix{1,2}=[7.70 0.00 0.00`.
- Lines 28: computes `inter.coupling.matrix{1,3}` using `inter.coupling.matrix{1,3}=[1.00 0.00 0.00`.
- Lines 33: computes `inter.relaxation` using `inter.relaxation={'t1_t2'}`.
- Lines 34: computes `inter.r1_rates` using `inter.r1_rates={1e3 1e4 1e4}`.
- Lines 35: computes `inter.r2_rates` using `inter.r2_rates={1e1 1e2 1e2}`.
- Lines 36: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 37: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 40: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 41: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 44: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 48: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 49: computes `parameters.electrons` using `parameters.electrons=1`.
- Lines 50: computes `parameters.nuclei` using `parameters.nuclei=[2 3]`.

## Implementation structure

- Mims ENDOR pulse sequence on BDPA with ideal electron pulses,
- reproducing Figure 10 from
- Calculation time: hours, much faster on GPU
- Isotopes
- Magnet field
- Interactions
- Relaxation theory
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `spin()`, `powder()`, `kfigure()`, `kylabel()`, `kxlabel()`.
