# examples/nmr_solids/cp_crystal_static_nh.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/cp_crystal_static_nh.m`
- Signature: `cp_crystal_static_nh()`
- Total lines: 51

## Purpose

1H-15N cross-polarisation experiment in the doubly rotating frame. Static single crystal simulation. Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: System specification; implemented by `sys.magnet=9.394`.
- Lines 14-15: Interactions; implemented by `inter.zeeman.scalar={0.00 0.00}`.
- Lines 20-21: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 24-25: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 28-29: Experiment parameters; implemented by `parameters.irr_powers=[5e4*ones(1,100)`.
- Lines 41-42: Simulation; implemented by `fid=crystal(spin_system,@cp_contact_hard,parameters,'nmr')`.
- Lines 44-45: Plotting; implemented by `time_axis=[0 cumsum(parameters.time_steps)]`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=9.394`.
- Lines 12: computes `sys.isotopes` using `sys.isotopes={'15N','1H'}`.
- Lines 15: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.00 0.00}`.
- Lines 16: computes `inter.coordinates` using `inter.coordinates={[0.00 0.00 0.00]`.
- Lines 18: computes `inter.temperature` using `inter.temperature=298`.
- Lines 21: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 22: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 25: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 29: computes `parameters.irr_powers` using `parameters.irr_powers=[5e4*ones(1,100)`.
- Lines 31-32: computes `parameters.irr_opers` using `parameters.irr_opers={operator(spin_system,'Ly','1H') operator(spin_system,'Lx','15N')}`.
- Lines 33-34: computes `parameters.exc_opers` using `parameters.exc_opers={operator(spin_system,'Lx','1H') operator(spin_system,'Ly','15N')}`.
- Lines 35: computes `parameters.coil` using `parameters.coil=state(spin_system,'Lx','15N')`.
- Lines 36: computes `parameters.time_steps` using `parameters.time_steps=1e-5*ones(1,100)`.
- Lines 37: computes `parameters.needs` using `parameters.needs={'aniso_eq'}`.
- Lines 38: computes `parameters.spins` using `parameters.spins={'15N'}`.
- Lines 39: computes `parameters.orientation` using `parameters.orientation=[pi/3 pi/4 pi/5]`.
- Lines 42: computes `fid` using `fid=crystal(spin_system,@cp_contact_hard,parameters,'nmr')`.
- Lines 45: computes `time_axis` using `time_axis=[0 cumsum(parameters.time_steps)]`.

## Implementation structure

- 1H-15N cross-polarisation experiment in the doubly rotating
- frame. Static single crystal simulation.
- Calculation time: seconds
- System specification
- Interactions
- Basis set
- Spinach housekeeping
- Experiment parameters
- Simulation
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `operator()`, `state()`, `crystal()`, `cumsum()`, `kfigure()`, `kylabel()`, `kxlabel()`.
