# examples/microfluidics/reacting_nmr.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/microfluidics/reacting_nmr.m`
- Signature: `reacting_nmr()`
- Total lines: 213

## Purpose

Non-linear reaction kinetics in combination with spin evolution (repeated pulse-acquire NMR) and relaxation (Redfield theory). Calculation time: hours, much faster on GPU.

## Physical / mathematical content

- Microfluidics examples. The coupled model is spin dynamics plus advection-diffusion-reaction transport on a mesh or regular grid. Numerical issues include finite-difference operators, mesh interpolation, and coupled reaction-flow evolution.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Import Diels-Alder cycloaddition; implemented by `[sys,inter,bas,kin]=dac_reaction()`.
- Lines 16-17: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 19-20: Greedy parallelisation; implemented by `sys.enable={'greedy'}`.
- Lines 22-23: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 26-27: Rate constants, mol/(L*s); implemented by `k1=0.5`.
- Lines 30-31: Cycloaddition reaction generator, including solvent; implemented by `K=@(t,x)(1i*[-k1*x(2)-k2*x(2) 0 0 0 0`.
- Lines 37-38: Kinetic time grid, 20 seconds; implemented by `chem_nsteps=200; chem_tmax=20`.
- Lines 42-43: Preallocate concentration trajectory; implemented by `chem_traj=zeros(5,chem_nsteps+1)`.
- Lines 45-46: Initial concentrations, mol/L; implemented by `chem_traj(:,1)=[0.6; 0.5; 0.0; 0.0; 18.1]`.
- Lines 48-49: Stage 1: concentration dynamics; implemented by `for n=1:chem_nsteps`.
- Lines 54-55: Plot concentrations, excluding solvent; implemented by `kfigure(); plot(chem_time_grid,real(chem_traj(1:4,:)))`.
- Lines 63-64: Interpolate concentrations as functions of time; implemented by `A=griddedInterpolant(chem_time_grid,chem_traj(1,:),'makima','none')`.
- Lines 69-70: Build chemical reaction generators; implemented by `G1=react_gen(spin_system,kin{1})`.
- Lines 73-77: Get concentration-weighted initial condition, no solvent; implemented by `eta= A(0)*state(spin_system,'Lz',spin_system.chem.parts{1}) +B(0)*state(spin_system,'Lz',spin_system.chem.parts{2}) +C(0)*state(spin_system,'Lz',spin_system.chem.parts{3…`.
- Lines 81-82: Preallocate the trajectory and get it started; implemented by `chem_traj=zeros([numel(eta) chem_nsteps+1]); chem_traj(:,1)=eta`.
- Lines 84-85: Run chemistry; implemented by `for n=1:chem_nsteps`.
- Lines 87-89: Keep the user informed; implemented by `report(spin_system,['chemistry time step ' int2str(n) '/' int2str(chem_nsteps)])`.
- Lines 91-92: Build the left interval edge composite evolution generator; implemented by `F_L=1i*k1*G1{1}*B(chem_time_grid(n))`.

### Control flow inferred from the code

- Line 49: `for` loop over `n=1:chem_nsteps`.
- Line 85: `for` loop over `n=1:chem_nsteps`.
- Line 132: `parfor` loop over `n=0:18`.
- Line 145: conditional branch on `ismember('gpu',spin_system.sys.enable)`.
- Line 160: `for` loop over `k=1:parameters.nsteps`.

### Key state/data transformations

- Lines 14: computes `[sys,inter,bas,kin]` using `[sys,inter,bas,kin]=dac_reaction()`.
- Lines 17: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 20: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 23: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 27: computes `k1` using `k1=0.5`.
- Lines 28: computes `k2` using `k2=0.1`.
- Lines 31: computes `K` using `K=@(t,x)(1i*[-k1*x(2)-k2*x(2) 0 0 0 0`.
- Lines 38: computes `chem_nsteps` using `chem_nsteps=200; chem_tmax=20`.
- Lines 39: computes `chem_dt` using `chem_dt=chem_tmax/chem_nsteps`.
- Lines 40: computes `chem_time_grid` using `chem_time_grid=linspace(0,chem_tmax,chem_nsteps+1)`.
- Lines 43: computes `chem_traj` using `chem_traj=zeros(5,chem_nsteps+1)`.
- Lines 46: computes `chem_traj(:,1)` using `chem_traj(:,1)=[0.6; 0.5; 0.0; 0.0; 18.1]`.
- Lines 50-51: computes `chem_traj(:,n+1)` using `chem_traj(:,n+1)=step(spin_system,{K,(n-1)*chem_dt,'LG4'}, chem_traj(:,n),chem_dt)`.
- Lines 64: computes `A` using `A=griddedInterpolant(chem_time_grid,chem_traj(1,:),'makima','none')`.
- Lines 65: computes `B` using `B=griddedInterpolant(chem_time_grid,chem_traj(2,:),'makima','none')`.
- Lines 66: computes `C` using `C=griddedInterpolant(chem_time_grid,chem_traj(3,:),'makima','none')`.
- Lines 67: computes `D` using `D=griddedInterpolant(chem_time_grid,chem_traj(4,:),'makima','none')`.
- Lines 70: computes `G1` using `G1=react_gen(spin_system,kin{1})`.

## Implementation structure

- Non-linear reaction kinetics in combination with spin evolution
- (repeated pulse-acquire NMR) and relaxation (Redfield theory).
- Calculation time: hours, much faster on GPU.
- Import Diels-Alder cycloaddition
- Magnet field
- Greedy parallelisation
- Spinach housekeeping
- Rate constants, mol/(L*s)
- Cycloaddition reaction generator, including solvent
- Kinetic time grid, 20 seconds
- Preallocate concentration trajectory
- Initial concentrations, mol/L

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `dac_reaction()`, `create()`, `basis()`, `chem_traj()`, `step()`, `kfigure()`, `kxlabel()`, `kylabel()`, `klegend()`, `griddedInterpolant()`, `react_gen()`, `state()`, `levelpop()`, `report()`, `int2str()`, `chem_time_grid()`.
