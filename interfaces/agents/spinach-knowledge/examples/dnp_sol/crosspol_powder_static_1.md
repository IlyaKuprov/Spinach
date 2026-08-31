# examples/dnp_sol/crosspol_powder_static_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/crosspol_powder_static_1.m`
- Signature: `crosspol_powder_static_1()`
- Total lines: 53

## Purpose

E-15N cross-polarization experiment in the doubly rotating frame. Static powder simulation. Calculation time: seconds

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: System specification; implemented by `sys.magnet=9.394`.
- Lines 14-15: Interactions; implemented by `inter.zeeman.scalar={0.00 2.0023193043622}`.
- Lines 20-21: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 24-25: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 28-29: Experiment parameters; implemented by `parameters.irr_powers=[5e4*ones(1,100)`.
- Lines 41-42: Simulation; implemented by `fid=powder(spin_system,@cp_contact_hard,parameters,'nmr')`.
- Lines 44-45: Time axis generation; implemented by `t_axis=[0 cumsum(parameters.time_steps)]`.
- Lines 47-48: Plotting; implemented by `kfigure(); plot(t_axis,real(fid)); kgrid`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=9.394`.
- Lines 12: computes `sys.isotopes` using `sys.isotopes={'15N','E'}`.
- Lines 15: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.00 2.0023193043622}`.
- Lines 16: computes `inter.coordinates` using `inter.coordinates={[0.00 0.00 0.00 ]`.
- Lines 18: computes `inter.temperature` using `inter.temperature=298`.
- Lines 21: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 22: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 25: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 29: computes `parameters.irr_powers` using `parameters.irr_powers=[5e4*ones(1,100)`.
- Lines 31-32: computes `parameters.irr_opers` using `parameters.irr_opers={operator(spin_system,'Ly','E'), operator(spin_system,'Lx','15N')}`.
- Lines 33-34: computes `parameters.exc_opers` using `parameters.exc_opers={operator(spin_system,'Lx','E'), operator(spin_system,'Ly','15N')}`.
- Lines 35: computes `parameters.needs` using `parameters.needs={'iso_eq'}`.
- Lines 36: computes `parameters.coil` using `parameters.coil=state(spin_system,'Lx','15N')`.
- Lines 37: computes `parameters.grid` using `parameters.grid='rep_2ang_6400pts_sph'`.
- Lines 38: computes `parameters.time_steps` using `parameters.time_steps=1e-5*ones(1,100)`.
- Lines 39: computes `parameters.spins` using `parameters.spins={'15N'}`.
- Lines 42: computes `fid` using `fid=powder(spin_system,@cp_contact_hard,parameters,'nmr')`.
- Lines 45: computes `t_axis` using `t_axis=[0 cumsum(parameters.time_steps)]`.

## Implementation structure

- E-15N cross-polarization experiment in the doubly rotating
- frame. Static powder simulation.
- Calculation time: seconds
- System specification
- Interactions
- Basis set
- Spinach housekeeping
- Experiment parameters
- Simulation
- Time axis generation
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `operator()`, `state()`, `powder()`, `cumsum()`, `kfigure()`, `kylabel()`, `kxlabel()`.
