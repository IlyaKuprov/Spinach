# examples/microfluidics/reacting_flow_nmr.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/microfluidics/reacting_flow_nmr.m`
- Signature: `reacting_flow_nmr()`
- Total lines: 252

## Purpose

Complete microfluidic simulation: diffusion, flow, two second- order chemical reactions, and NMR detection in a narrow strip of the chip where the coil is assumed to be located. Calculation time: days, much faster on GPU.

## Physical / mathematical content

- Microfluidics examples. The coupled model is spin dynamics plus advection-diffusion-reaction transport on a mesh or regular grid. Numerical issues include finite-difference operators, mesh interpolation, and coupled reaction-flow evolution.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Import Diels-Alder cycloaddition; implemented by `[sys,inter,bas,kin]=dac_reaction()`.
- Lines 19-20: Import hydrodynamics information; implemented by `comsol.mesh_file='chip_mesh.txt'`.
- Lines 31-32: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 34-35: This needs a GPU; implemented by `sys.enable={'greedy'}`.
- Lines 37-38: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 42-45: % Concentration dynamics stage; implemented by `k1=2.0`.
- Lines 44-45: Rate constants, mol/(L*s); implemented by `k1=2.0`.
- Lines 48-49: Cycloaddition reaction generator, including solvent; implemented by `K=@(x)([-k1*x(2)-k2*x(2) 0 0 0 0`.
- Lines 55-56: Strong diffusion; implemented by `parameters.diff=1e-7`.
- Lines 58-59: Timing parameters; implemented by `chem_dt=20; chem_nsteps=501`.
- Lines 61-62: Get diffusion and flow generator; implemented by `GF=flow_gen(spin_system,parameters)`.
- Lines 64-65: Concentration trajectory preallocation and the initial state; implemented by `chem_traj=zeros(5,spin_system.mesh.vor.ncells,chem_nsteps+1)`.
- Lines 68-69: Time evolution loop; implemented by `for n=1:chem_nsteps`.
- Lines 71-73: Keep the user informed; implemented by `report(spin_system,['chemistry + hydrodynamics time step ' int2str(n) '/' int2str(chem_nsteps)])`.
- Lines 75-76: Build a kinetics generator in each cell; implemented by `GK=zeros([5 5 spin_system.mesh.vor.ncells],'like',1i)`.
- Lines 81-82: Assemble the evolution generator; implemented by `G=1i*sp_block_diag(GK)+1i*kron(GF,speye(5))`.
- Lines 84-85: Take the time step; implemented by `c_curr=chem_traj(:,:,n); c_curr=c_curr(:)`.
- Lines 91-92: Get chemistry time grid; implemented by `chem_time_grid=linspace(0,chem_dt*chem_nsteps,chem_nsteps+1)`.

### Control flow inferred from the code

- Line 69: `for` loop over `n=1:chem_nsteps`.
- Line 77: `parfor` loop over `k=1:spin_system.mesh.vor.ncells`.
- Line 97: `parfor` loop over `n=1:spin_system.mesh.vor.ncells`.
- Line 143: conditional branch on `ismember('gpu',sys.enable)`.
- Line 160: `parfor` loop over `j=1:numel(n_vals)`.
- Line 168: `for` loop over `k=1:spin_system.mesh.vor.ncells`.
- Line 187: `for` loop over `k=1:parameters.nsteps`.
- Line 196: `for` loop over `m=1:spin_system.mesh.vor.ncells`.

### Key state/data transformations

- Lines 17: computes `[sys,inter,bas,kin]` using `[sys,inter,bas,kin]=dac_reaction()`.
- Lines 20: computes `comsol.mesh_file` using `comsol.mesh_file='chip_mesh.txt'`.
- Lines 21: computes `comsol.velo_file` using `comsol.velo_file='chip_velo.txt'`.
- Lines 22: computes `comsol.crop` using `comsol.crop={[286.8 287.5],[576.0 579.0]}`.
- Lines 23-28: computes `comsol.inactivate` using `comsol.inactivate=[9 10 19 30 20 25 14 13 3372 3373 3380 3381 3382 3386 3169 3185 3201 3054 3077 3055 3053 3078 3186 3168 875 899 897 877 876 860 858 885 859 883]`.
- Lines 29: computes `mesh` using `mesh=comsol_import(comsol)`.
- Lines 32: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 35: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 38: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 40: computes `spin_system.mesh` using `spin_system.mesh=mesh`.
- Lines 45: computes `k1` using `k1=2.0`.
- Lines 46: computes `k2` using `k2=1.0`.
- Lines 49: computes `K` using `K=@(x)([-k1*x(2)-k2*x(2) 0 0 0 0`.
- Lines 56: computes `parameters.diff` using `parameters.diff=1e-7`.
- Lines 59: computes `chem_dt` using `chem_dt=20; chem_nsteps=501`.
- Lines 62: computes `GF` using `GF=flow_gen(spin_system,parameters)`.
- Lines 65: computes `chem_traj` using `chem_traj=zeros(5,spin_system.mesh.vor.ncells,chem_nsteps+1)`.
- Lines 66: computes `chem_traj(1,1240,1)` using `chem_traj(1,1240,1)=0.50; chem_traj(2,1246,1)=0.25`.

## Implementation structure

- Complete microfluidic simulation: diffusion, flow, two second-
- order chemical reactions, and NMR detection in a narrow strip
- of the chip where the coil is assumed to be located.
- Calculation time: days, much faster on GPU.
- Import Diels-Alder cycloaddition
- Import hydrodynamics information
- Magnet field
- This needs a GPU
- Spinach housekeeping
- % Concentration dynamics stage
- Rate constants, mol/(L*s)
- Cycloaddition reaction generator, including solvent

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `dac_reaction()`, `comsol_import()`, `create()`, `basis()`, `flow_gen()`, `chem_traj()`, `report()`, `int2str()`, `sp_block_diag()`, `speye()`, `c_curr()`, `step()`, `griddedInterpolant()`, `squeeze()`, `react_gen()`, `double()`.
