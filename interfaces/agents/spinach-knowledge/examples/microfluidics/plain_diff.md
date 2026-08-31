# examples/microfluidics/plain_diff.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/microfluidics/plain_diff.m`
- Signature: `plain_diff()`
- Total lines: 115

## Purpose

Simple diffusion simulation without spin dynamics. Longitudinal magnetisation is tracked as a function of time.

## Physical / mathematical content

- Microfluidics examples. The coupled model is spin dynamics plus advection-diffusion-reaction transport on a mesh or regular grid. Numerical issues include finite-difference operators, mesh interpolation, and coupled reaction-flow evolution.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Import hydrodynamics information; implemented by `comsol.mesh_file='chip_mesh.txt'`.
- Lines 23-24: One proton; implemented by `sys.magnet=14.1`.
- Lines 27-28: Chemical shift (water); implemented by `inter.zeeman.scalar={0.0}`.
- Lines 30-31: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 34-35: Algorithmic switches; implemented by `sys.disable={'trajlevel'}`.
- Lines 37-38: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 42-43: Initial condition: Lz in one cell in the middle; implemented by `parameters.rho0_ph{1}=zeros(spin_system.mesh.vor.ncells,1)`.
- Lines 47-48: Detection state: Lz in all cells; implemented by `parameters.coil_ph{1}=ones(spin_system.mesh.vor.ncells,1)`.
- Lines 51-52: Sequence and timing parameters; implemented by `parameters.spins={'1H'}`.
- Lines 57-58: Set assumptions; implemented by `spin_system=assume(spin_system,'nmr')`.
- Lines 60-61: Same Hamiltonian everywhere; implemented by `H=hamiltonian(spin_system)`.
- Lines 66-67: Same relaxation everywhere; implemented by `parameters.R_op={relaxation(spin_system)}`.
- Lines 70-71: Just diffusion; implemented by `spin_system.mesh.u=0*spin_system.mesh.u`.
- Lines 75-76: Drainage in the distal pipe; implemented by `drainage=zeros(2659,1)`.
- Lines 81-82: Get the trajectory of simple flow; implemented by `traj=meshflow(spin_system,@simple_flow,parameters)`.
- Lines 84-85: Extract the observable quantity; implemented by `coil=state(spin_system,'Lz','1H')`.
- Lines 88-89: Make a figure; implemented by `kfigure(); scale_figure([1.5 1.5])`.
- Lines 92-93: Set Z axis extents; implemented by `spin_system.mesh.zext=[-0.02 0.02]`.

### Control flow inferred from the code

- Line 96: `for` loop over `n=1:size(traj,2)`.

### Key state/data transformations

- Lines 12: computes `comsol.mesh_file` using `comsol.mesh_file='chip_mesh.txt'`.
- Lines 13: computes `comsol.velo_file` using `comsol.velo_file='chip_velo.txt'`.
- Lines 14: computes `comsol.crop` using `comsol.crop={[286.8 287.5],[576.0 579.0]}`.
- Lines 15-20: computes `comsol.inactivate` using `comsol.inactivate=[9 10 19 30 20 25 14 13 3372 3373 3380 3381 3382 3386 3169 3185 3201 3054 3077 3055 3053 3078 3186 3168 875 899 897 877 876 860 858 885 859 883]`.
- Lines 21: computes `mesh` using `mesh=comsol_import(comsol)`.
- Lines 24: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 25: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 28: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0}`.
- Lines 31: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 32: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 35: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 38: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 40: computes `spin_system.mesh` using `spin_system.mesh=mesh`.
- Lines 43: computes `parameters.rho0_ph{1}` using `parameters.rho0_ph{1}=zeros(spin_system.mesh.vor.ncells,1)`.
- Lines 44: computes `parameters.rho0_ph{1}([1240 1246])` using `parameters.rho0_ph{1}([1240 1246])=0.5`.
- Lines 45: computes `parameters.rho0_st{1}` using `parameters.rho0_st{1}=state(spin_system,'Lz','1H')`.
- Lines 48: computes `parameters.coil_ph{1}` using `parameters.coil_ph{1}=ones(spin_system.mesh.vor.ncells,1)`.
- Lines 49: computes `parameters.coil_st{1}` using `parameters.coil_st{1}=state(spin_system,'Lz','1H')`.

## Implementation structure

- Simple diffusion simulation without spin dynamics. Longitudinal
- magnetisation is tracked as a function of time.
- Import hydrodynamics information
- One proton
- Chemical shift (water)
- Basis set
- Algorithmic switches
- Spinach housekeeping
- Initial condition: Lz in one cell in the middle
- Detection state: Lz in all cells
- Sequence and timing parameters
- Set assumptions

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `comsol_import()`, `create()`, `basis()`, `state()`, `assume()`, `hamiltonian()`, `frqoffset()`, `relaxation()`, `drainage()`, `speye()`, `meshflow()`, `fpl2phan()`, `traj()`, `kfigure()`, `scale_figure()`, `camproj()`.
