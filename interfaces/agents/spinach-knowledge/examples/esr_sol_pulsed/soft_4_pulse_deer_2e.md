# examples/esr_sol_pulsed/soft_4_pulse_deer_2e.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/soft_4_pulse_deer_2e.m`
- Signature: `soft_4_pulse_deer_2e()`
- Total lines: 80

## Purpose

Four-pulse DEER simulation for a two-electron system. Soft pulses are simulated using the Fokker-Planck formalism. Calculation time: minutes

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Magnet field; implemented by `sys.magnet=0.3451805`.
- Lines 13-14: Isotopes; implemented by `sys.isotopes={'E','E'}`.
- Lines 16-17: Zeeman interactions; implemented by `inter.zeeman.eigs=cell(1,2)`.
- Lines 24-26: Spin-orbit corrections to the DD couplings; implemented by `sys.enable={'sodd'}`.
- Lines 28-29: Coordinates (Angstrom); implemented by `inter.coordinates=cell(2,1)`.
- Lines 33-34: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 37-38: Algorithmic options; implemented by `sys.disable={'trajlevel'}`.
- Lines 40-41: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 44-45: Sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 52-53: EPR parameters; implemented by `parameters.offset=-4e8`.
- Lines 62-63: DEER pulse parameters; implemented by `parameters.pulse_rnk=[2 2 2 2]`.
- Lines 69-70: DEER echo timing parameters; implemented by `parameters.p1_p2_gap=0.5e-6`.
- Lines 76-77: Simulation and plotting; implemented by `deer_4p_soft_diag(spin_system,parameters)`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=0.3451805`.
- Lines 14: computes `sys.isotopes` using `sys.isotopes={'E','E'}`.
- Lines 17: computes `inter.zeeman.eigs` using `inter.zeeman.eigs=cell(1,2)`.
- Lines 18: computes `inter.zeeman.euler` using `inter.zeeman.euler=cell(1,2)`.
- Lines 19: computes `inter.zeeman.eigs{1}` using `inter.zeeman.eigs{1}=[2.284 2.123 2.075]`.
- Lines 20: computes `inter.zeeman.euler{1}` using `inter.zeeman.euler{1}=[135 90 45]*(pi/180)`.
- Lines 21: computes `inter.zeeman.eigs{2}` using `inter.zeeman.eigs{2}=[2.035 2.013 1.975]`.
- Lines 22: computes `inter.zeeman.euler{2}` using `inter.zeeman.euler{2}=[30 60 120]*(pi/180)`.
- Lines 26: computes `sys.enable` using `sys.enable={'sodd'}`.
- Lines 29: computes `inter.coordinates` using `inter.coordinates=cell(2,1)`.
- Lines 30: computes `inter.coordinates{1}` using `inter.coordinates{1}=[0.00 0.00 0.00]`.
- Lines 31: computes `inter.coordinates{2}` using `inter.coordinates{2}=[20.00 0.00 0.00]`.
- Lines 34: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 35: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 38: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 41: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 45: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 46: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lz','E')`.

## Implementation structure

- Four-pulse DEER simulation for a two-electron system. Soft
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

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `deer_4p_soft_diag()`.
