# examples/spin_chemistry/singlet_yield_anisotropy_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/spin_chemistry/singlet_yield_anisotropy_2.m`
- Signature: `singlet_yield_anisotropy_2()`
- Total lines: 73

## Purpose

Singlet yield anisotropy calculation for a radical pair using exponential recombination kinetics model. Calculation time: seconds

## Physical / mathematical content

- Spin-chemistry examples. These scripts treat radical pairs, recombination channels, chemically induced dynamic nuclear polarisation, and magnetic-field effects. The theory combines spin-selective kinetics with singlet-triplet interconversion.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Unit magnet (field sweep); implemented by `sys.magnet=1`.
- Lines 15-16: Isotopes; implemented by `sys.isotopes={'E','E','14N','14N','1H'}`.
- Lines 18-19: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 22-23: Rotation matrices; implemented by `R1=[0.4380 0.8655 -0.2432; 0.8981 -0.4097 0.1595; -0.0384 0.2883 0.9568]`.
- Lines 27-28: Interaction eigenvalues; implemented by `A1=[-1.049 0 0; 0 -0.996 0; 0 0 13.826]`.
- Lines 32-33: Coupling tensors; implemented by `inter.coupling.matrix=cell(5)`.
- Lines 38-39: Zeeman interactions; implemented by `inter.zeeman.scalar={2.0023 2.0023 0 0 0}`.
- Lines 41-42: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 45-46: Sequence parameters; implemented by `parameters.npoints=1`.
- Lines 55-56: Simulation; implemented by `[yield,grid]=powder(spin_system,@rydmr_exp,parameters,'labframe')`.
- Lines 58-59: Preprocessing; implemented by `yield=cell2mat(yield)`.
- Lines 62-63: Plotting; implemented by `hull=get_hull(grid.betas,grid.gammas)`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=1`.
- Lines 16: computes `sys.isotopes` using `sys.isotopes={'E','E','14N','14N','1H'}`.
- Lines 19: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 20: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 23: computes `R1` using `R1=[0.4380 0.8655 -0.2432; 0.8981 -0.4097 0.1595; -0.0384 0.2883 0.9568]`.
- Lines 24: computes `R2` using `R2=[0.9703 -0.2207 0.0992; 0.2383 0.9426 -0.2340; -0.0419 0.2506 0.9672]`.
- Lines 25: computes `R3` using `R3=[0.9819 0.1883 -0.0203; -0.0348 0.2850 0.9579; -0.1861 0.9398 -0.2864]`.
- Lines 28: computes `A1` using `A1=[-1.049 0 0; 0 -0.996 0; 0 0 13.826]`.
- Lines 29: computes `A2` using `A2=[-0.305 0 0; 0 -0.222 0; 0 0 6.872]`.
- Lines 30: computes `A3` using `A3=[-13.850 0 0; 0 -9.372 0; 0 0 0.143]`.
- Lines 33: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(5)`.
- Lines 34: computes `inter.coupling.matrix{1,3}` using `inter.coupling.matrix{1,3}=1e6*gauss2mhz(R1*A1*R1')`.
- Lines 35: computes `inter.coupling.matrix{2,4}` using `inter.coupling.matrix{2,4}=1e6*gauss2mhz(R2*A2*R2')`.
- Lines 36: computes `inter.coupling.matrix{1,5}` using `inter.coupling.matrix{1,5}=1e6*gauss2mhz(R3*A3*R3')`.
- Lines 39: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={2.0023 2.0023 0 0 0}`.
- Lines 42: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 46: computes `parameters.npoints` using `parameters.npoints=1`.
- Lines 47: computes `parameters.fields` using `parameters.fields=50e-6`.

## Implementation structure

- Singlet yield anisotropy calculation for a radical pair
- using exponential recombination kinetics model.
- Calculation time: seconds
- Unit magnet (field sweep)
- Isotopes
- Basis set
- Rotation matrices
- Interaction eigenvalues
- Coupling tensors
- Zeeman interactions
- Spinach housekeeping
- Sequence parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gauss2mhz()`, `create()`, `basis()`, `powder()`, `cell2mat()`, `get_hull()`, `kfigure()`, `trisurf()`.
