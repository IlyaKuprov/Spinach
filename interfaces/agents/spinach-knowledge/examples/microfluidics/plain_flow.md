# examples/microfluidics/plain_flow.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/microfluidics/plain_flow.m`
- Signature: `plain_flow()`
- Total lines: 116

## Purpose

Simple flow simulation with no dynamics in the spin subspace: longitudinal magnetisation is tracked as a function of time af- ter injection into the flow field imported from COMSOL with a diffusion term also present. The tail of the pipe has drainage terms set up using a kinetics superoperator phantom.

## Physical / mathematical content

- Microfluidics examples. The coupled model is spin dynamics plus advection-diffusion-reaction transport on a mesh or regular grid. Numerical issues include finite-difference operators, mesh interpolation, and coupled reaction-flow evolution.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Import hydrodynamics information; implemented by `comsol.mesh_file='chip_mesh.txt'`.
- Lines 26-27: One proton; implemented by `sys.magnet=14.1`.
- Lines 30-31: Chemical shift (water); implemented by `inter.zeeman.scalar={0.0}`.
- Lines 33-34: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 37-38: Algorithmic switches; implemented by `sys.disable={'trajlevel'}`.
- Lines 40-41: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 45-46: Initial condition: Lz in a few cells; implemented by `parameters.rho0_ph{1}=zeros(spin_system.mesh.vor.ncells,1)`.
- Lines 50-51: Detection state: Lz in all cells; implemented by `parameters.coil_ph{1}=ones(spin_system.mesh.vor.ncells,1)`.
- Lines 54-55: Sequence and timing parameters; implemented by `parameters.spins={'1H'}`.
- Lines 60-61: Set assumptions; implemented by `spin_system=assume(spin_system,'nmr')`.
- Lines 63-64: Same Hamiltonian everywhere; implemented by `H=hamiltonian(spin_system)`.
- Lines 69-70: Same relaxation everywhere; implemented by `parameters.R_op={relaxation(spin_system)}`.
- Lines 73-74: Diffusion coefficient; implemented by `parameters.diff=1e-7`.
- Lines 76-77: Drainage in the distal pipe; implemented by `drainage=zeros(2659,1)`.
- Lines 82-83: Get the trajectory of simple flow; implemented by `traj=meshflow(spin_system,@simple_flow,parameters)`.
- Lines 85-86: Extract the observable quantity; implemented by `coil=state(spin_system,'Lz','1H')`.
- Lines 89-90: Make a figure; implemented by `kfigure(); scale_figure([1.5 1.5])`.
- Lines 93-94: Set Z axis extents; implemented by `spin_system.mesh.zext=[-0.02 0.02]`.

### Control flow inferred from the code

- Line 97: `for` loop over `n=1:size(traj,2)`.

### Key state/data transformations

- Lines 15: computes `comsol.mesh_file` using `comsol.mesh_file='chip_mesh.txt'`.
- Lines 16: computes `comsol.velo_file` using `comsol.velo_file='chip_velo.txt'`.
- Lines 17: computes `comsol.crop` using `comsol.crop={[286.8 287.5],[576.0 579.0]}`.
- Lines 18-23: computes `comsol.inactivate` using `comsol.inactivate=[9 10 19 30 20 25 14 13 3372 3373 3380 3381 3382 3386 3169 3185 3201 3054 3077 3055 3053 3078 3186 3168 875 899 897 877 876 860 858 885 859 883]`.
- Lines 24: computes `mesh` using `mesh=comsol_import(comsol)`.
- Lines 27: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 28: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 31: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0}`.
- Lines 34: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 35: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 38: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 41: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 43: computes `spin_system.mesh` using `spin_system.mesh=mesh`.
- Lines 46: computes `parameters.rho0_ph{1}` using `parameters.rho0_ph{1}=zeros(spin_system.mesh.vor.ncells,1)`.
- Lines 47: computes `parameters.rho0_ph{1}(140:160)` using `parameters.rho0_ph{1}(140:160)=2`.
- Lines 48: computes `parameters.rho0_st{1}` using `parameters.rho0_st{1}=state(spin_system,'Lz','1H')`.
- Lines 51: computes `parameters.coil_ph{1}` using `parameters.coil_ph{1}=ones(spin_system.mesh.vor.ncells,1)`.
- Lines 52: computes `parameters.coil_st{1}` using `parameters.coil_st{1}=state(spin_system,'Lz','1H')`.

## Implementation structure

- Simple flow simulation with no dynamics in the spin subspace:
- longitudinal magnetisation is tracked as a function of time af-
- ter injection into the flow field imported from COMSOL with a
- diffusion term also present. The tail of the pipe has drainage
- terms set up using a kinetics superoperator phantom.
- Import hydrodynamics information
- One proton
- Chemical shift (water)
- Basis set
- Algorithmic switches
- Spinach housekeeping
- Initial condition: Lz in a few cells

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `comsol_import()`, `create()`, `basis()`, `state()`, `assume()`, `hamiltonian()`, `frqoffset()`, `relaxation()`, `drainage()`, `speye()`, `meshflow()`, `fpl2phan()`, `traj()`, `kfigure()`, `scale_figure()`, `camproj()`.
