# examples/esr_sol_pulsed/soft_3_pulse_deer_3e.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/soft_3_pulse_deer_3e.m`
- Signature: `soft_3_pulse_deer_3e()`
- Total lines: 82

## Purpose

Three-pulse DEER simulation for a three-electron system. Soft pulses are simulated using the Fokker-Planck formalism. Calculation time: minutes

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Magnet field; implemented by `sys.magnet=0.3451805`.
- Lines 13-14: Isotopes; implemented by `sys.isotopes={'E','E','E'}`.
- Lines 16-17: Zeeman interactions; implemented by `inter.zeeman.eigs=cell(1,2)`.
- Lines 26-28: Spin-orbit corrections to the DD couplings; implemented by `sys.enable={'sodd'}`.
- Lines 30-31: Coordinates (Angstrom); implemented by `inter.coordinates=cell(3,1)`.
- Lines 36-37: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 40-41: Algorithmic options; implemented by `sys.disable={'trajlevel'}`.
- Lines 43-44: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 47-48: Sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 55-56: EPR parameters; implemented by `parameters.offset=-4e8`.
- Lines 65-66: DEER pulse parameters; implemented by `parameters.pulse_rnk=[2 2 2]`.
- Lines 72-73: DEER echo timing parameters; implemented by `parameters.p1_p3_gap=1e-6`.
- Lines 78-79: Simulation and plotting; implemented by `deer_3p_soft_diag(spin_system,parameters)`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=0.3451805`.
- Lines 14: computes `sys.isotopes` using `sys.isotopes={'E','E','E'}`.
- Lines 17: computes `inter.zeeman.eigs` using `inter.zeeman.eigs=cell(1,2)`.
- Lines 18: computes `inter.zeeman.euler` using `inter.zeeman.euler=cell(1,2)`.
- Lines 19: computes `inter.zeeman.eigs{1}` using `inter.zeeman.eigs{1}=[2.284 2.123 2.075]`.
- Lines 20: computes `inter.zeeman.euler{1}` using `inter.zeeman.euler{1}=[135 90 45]*(pi/180)`.
- Lines 21: computes `inter.zeeman.eigs{2}` using `inter.zeeman.eigs{2}=[2.035 2.013 1.975]`.
- Lines 22: computes `inter.zeeman.euler{2}` using `inter.zeeman.euler{2}=[30 60 120]*(pi/180)`.
- Lines 23: computes `inter.zeeman.eigs{3}` using `inter.zeeman.eigs{3}=[1.935 1.895 1.895]`.
- Lines 24: computes `inter.zeeman.euler{3}` using `inter.zeeman.euler{3}=[60 40 20]*(pi/180)`.
- Lines 28: computes `sys.enable` using `sys.enable={'sodd'}`.
- Lines 31: computes `inter.coordinates` using `inter.coordinates=cell(3,1)`.
- Lines 32: computes `inter.coordinates{1}` using `inter.coordinates{1}=[0.00 0.00 0.00]`.
- Lines 33: computes `inter.coordinates{2}` using `inter.coordinates{2}=[20.00 0.00 0.00]`.
- Lines 34: computes `inter.coordinates{3}` using `inter.coordinates{3}=[0.00 0.00 20.00]`.
- Lines 37: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 38: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 41: computes `sys.disable` using `sys.disable={'trajlevel'}`.

## Implementation structure

- Three-pulse DEER simulation for a three-electron system. Soft
- pulses are simulated using the Fokker-Planck formalism.
- Calculation time: minutes
- Magnet field
- Isotopes
- Zeeman interactions
- Spin-orbit corrections
- to the DD couplings
- Coordinates (Angstrom)
- Basis set
- Algorithmic options
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `deer_3p_soft_diag()`.
