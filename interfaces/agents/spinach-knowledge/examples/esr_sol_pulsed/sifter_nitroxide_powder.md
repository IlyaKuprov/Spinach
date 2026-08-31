# examples/esr_sol_pulsed/sifter_nitroxide_powder.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/sifter_nitroxide_powder.m`
- Signature: `sifter_nitroxide_powder()`
- Total lines: 72

## Purpose

An example of the SIFTER sequence. Calculation time: minutes.

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Magnet field; implemented by `sys.magnet=0.33`.
- Lines 13-14: System specification; implemented by `sys.isotopes={'E','E','14N','14N'}`.
- Lines 16-17: Zeeman interactions; implemented by `inter.zeeman.eigs=cell(1,4)`.
- Lines 24-25: Coordinates for inter-electron DD; implemented by `inter.coordinates={[0 0 0]; [0 0 20]; []; [] }`.
- Lines 27-28: Hyperfine couplings; implemented by `inter.coupling.eigs=cell(4,4)`.
- Lines 35-36: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 40-41: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 44-45: Set the sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 55-56: Simulation and time axis generation; implemented by `fid=imag(powder(spin_system,@sifter,parameters,'esr'))`.
- Lines 62-63: Plotting; implemented by `kfigure(); scale_figure([1.50 0.75])`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=0.33`.
- Lines 14: computes `sys.isotopes` using `sys.isotopes={'E','E','14N','14N'}`.
- Lines 17: computes `inter.zeeman.eigs` using `inter.zeeman.eigs=cell(1,4)`.
- Lines 18: computes `inter.zeeman.euler` using `inter.zeeman.euler=cell(1,4)`.
- Lines 19: computes `inter.zeeman.eigs{1,1}` using `inter.zeeman.eigs{1,1}=[2.0087 2.0058 2.0018]`.
- Lines 20: computes `inter.zeeman.eigs{1,2}` using `inter.zeeman.eigs{1,2}=[2.0087 2.0058 2.0018]`.
- Lines 21: computes `inter.zeeman.euler{1,1}` using `inter.zeeman.euler{1,1}=[0 0 0]`.
- Lines 22: computes `inter.zeeman.euler{1,2}` using `inter.zeeman.euler{1,2}=[0 0 0]`.
- Lines 25: computes `inter.coordinates` using `inter.coordinates={[0 0 0]; [0 0 20]; []; [] }`.
- Lines 28: computes `inter.coupling.eigs` using `inter.coupling.eigs=cell(4,4)`.
- Lines 29: computes `inter.coupling.euler` using `inter.coupling.euler=cell(4,4)`.
- Lines 30: computes `inter.coupling.eigs{1,3}` using `inter.coupling.eigs{1,3}=[19.8977 20.1780 102.8516]*1e6`.
- Lines 31: computes `inter.coupling.eigs{2,4}` using `inter.coupling.eigs{2,4}=[19.8977 20.1780 102.8516]*1e6`.
- Lines 32: computes `inter.coupling.euler{1,3}` using `inter.coupling.euler{1,3}=[0 0 0]`.
- Lines 33: computes `inter.coupling.euler{2,4}` using `inter.coupling.euler{2,4}=[0 0 0]`.
- Lines 36: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 37: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 38: computes `bas.longitudinals` using `bas.longitudinals={'14N'}`.

## Implementation structure

- An example of the SIFTER sequence.
- Calculation time: minutes.
- Magnet field
- System specification
- Zeeman interactions
- Coordinates for inter-electron DD
- Hyperfine couplings
- Basis set
- Spinach housekeeping
- Set the sequence parameters
- Simulation and time axis generation
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `powder()`, `kfigure()`, `scale_figure()`, `subplot()`, `kxlabel()`, `kylabel()`, `ktitle()`.
