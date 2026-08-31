# examples/nmr_solids/cp_crystal_static_nhh.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/cp_crystal_static_nhh.m`
- Signature: `cp_crystal_static_nhh()`
- Total lines: 67

## Purpose

Cross-polarisation experiment in the doubly rotating frame. A single nitrogen-15 in a bath of 8 protons scattered on a 2 Angstrom radius sphere around it. Static single crystal simulation in a full Liouvil- le space (here necessary because this is not a powder and everything interacts with everything). Calculation time: minutes on a Tesla A100, much longer on CPU.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: System specification; implemented by `sys.magnet=9.394`.
- Lines 19-21: Interactions; implemented by `inter.zeeman.scalar={0.6745 1.0368 -0.1495 -0.3171 0.7233 0.4882 3.2662 0.0317 0.0000}`.
- Lines 33-34: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 37-38: This needs a GPU; implemented by `sys.enable={'greedy'}`.
- Lines 40-41: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 44-45: Experiment parameters; implemented by `parameters.irr_powers=[5e4*ones(1,100)`.
- Lines 57-58: Simulation; implemented by `fid=crystal(spin_system,@cp_contact_hard,parameters,'nmr')`.
- Lines 60-61: Plotting; implemented by `time_axis=[0 cumsum(parameters.time_steps)]`.

### Key state/data transformations

- Lines 15: computes `sys.magnet` using `sys.magnet=9.394`.
- Lines 16-17: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H', '1H','1H','1H','1H','15N'}`.
- Lines 20-21: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.6745 1.0368 -0.1495 -0.3171 0.7233 0.4882 3.2662 0.0317 0.0000}`.
- Lines 22: computes `inter.coordinates` using `inter.coordinates={[-2.51887819 -0.99807636 0.87365165]`.
- Lines 31: computes `inter.temperature` using `inter.temperature=298`.
- Lines 34: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 35: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 38: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 41: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 45: computes `parameters.irr_powers` using `parameters.irr_powers=[5e4*ones(1,100)`.
- Lines 47-48: computes `parameters.irr_opers` using `parameters.irr_opers={operator(spin_system,'Ly','1H') operator(spin_system,'Lx','15N')}`.
- Lines 49-50: computes `parameters.exc_opers` using `parameters.exc_opers={operator(spin_system,'Lx','1H') operator(spin_system,'Ly','15N')}`.
- Lines 51: computes `parameters.coil` using `parameters.coil=state(spin_system,'Lx','15N')`.
- Lines 52: computes `parameters.time_steps` using `parameters.time_steps=1e-5*ones(1,100)`.
- Lines 53: computes `parameters.needs` using `parameters.needs={'aniso_eq'}`.
- Lines 54: computes `parameters.spins` using `parameters.spins={'15N'}`.
- Lines 55: computes `parameters.orientation` using `parameters.orientation=[pi/3 pi/4 pi/5]`.
- Lines 58: computes `fid` using `fid=crystal(spin_system,@cp_contact_hard,parameters,'nmr')`.

## Implementation structure

- Cross-polarisation experiment in the doubly rotating frame. A single
- nitrogen-15 in a bath of 8 protons scattered on a 2 Angstrom radius
- sphere around it. Static single crystal simulation in a full Liouvil-
- le space (here necessary because this is not a powder and everything
- interacts with everything).
- Calculation time: minutes on a Tesla A100, much longer on CPU.
- System specification
- Interactions
- Basis set
- This needs a GPU
- Spinach housekeeping
- Experiment parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `operator()`, `state()`, `crystal()`, `cumsum()`, `kfigure()`, `kylabel()`, `kxlabel()`.
