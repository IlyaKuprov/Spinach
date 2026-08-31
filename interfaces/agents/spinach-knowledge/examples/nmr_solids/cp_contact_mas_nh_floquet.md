# examples/nmr_solids/cp_contact_mas_nh_floquet.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/cp_contact_mas_nh_floquet.m`
- Signature: `cp_contact_mas_nh_floquet()`
- Total lines: 60

## Purpose

Cross-polarisation experiment in the doubly rotating frame. A single nitrogen-15 and a single proton. Spinning powder simulation starting from the thermal equilibrium using Floquet formalism. Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file relies on Floquet theory, where periodic time dependence is lifted into an enlarged block representation that converts time-periodic dynamics into a time-independent eigenproblem.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: System specification; implemented by `sys.magnet=9.394`.
- Lines 16-17: Interactions; implemented by `inter.zeeman.scalar={0.00 0.00}`.
- Lines 22-23: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 26-27: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 30-31: Relevant operators; implemented by `Nx=operator(spin_system,'Lx','15N')`.
- Lines 36-37: MAS parameters; implemented by `parameters.rate=10000`.
- Lines 50-51: Simulation; implemented by `fid=floquet(spin_system,@cp_contact_hard,parameters,'nmr')`.
- Lines 53-54: Plot the answer; implemented by `time_axis=[0 cumsum(parameters.time_steps)]`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=9.394`.
- Lines 14: computes `sys.isotopes` using `sys.isotopes={'15N','1H'}`.
- Lines 17: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.00 0.00}`.
- Lines 18: computes `inter.coordinates` using `inter.coordinates={[0.00 0.00 0.00]`.
- Lines 20: computes `inter.temperature` using `inter.temperature=298`.
- Lines 23: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 24: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 27: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 31: computes `Nx` using `Nx=operator(spin_system,'Lx','15N')`.
- Lines 32: computes `Ny` using `Ny=operator(spin_system,'Ly','15N')`.
- Lines 33: computes `Hx` using `Hx=operator(spin_system,'Lx','1H')`.
- Lines 34: computes `Hy` using `Hy=operator(spin_system,'Ly','1H')`.
- Lines 37: computes `parameters.rate` using `parameters.rate=10000`.
- Lines 38: computes `parameters.axis` using `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 39: computes `parameters.max_rank` using `parameters.max_rank=4`.
- Lines 40: computes `parameters.spins` using `parameters.spins={'1H','15N'}`.
- Lines 41: computes `parameters.irr_powers` using `parameters.irr_powers=[5e4*ones(1,100)`.
- Lines 43: computes `parameters.irr_opers` using `parameters.irr_opers={Hy Nx}`.

## Implementation structure

- Cross-polarisation experiment in the doubly rotating frame. A single
- nitrogen-15 and a single proton. Spinning powder simulation starting
- from the thermal equilibrium using Floquet formalism.
- Calculation time: seconds
- System specification
- Interactions
- Basis set
- Spinach housekeeping
- Relevant operators
- MAS parameters
- Simulation
- Plot the answer

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `operator()`, `state()`, `floquet()`, `cumsum()`, `kfigure()`, `kylabel()`, `kxlabel()`.
