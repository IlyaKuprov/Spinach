# examples/spin_chemistry/singlet_yield_anisotropy_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/spin_chemistry/singlet_yield_anisotropy_1.m`
- Signature: `singlet_yield_anisotropy_1()`
- Total lines: 57

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
- Lines 15-16: System specification; implemented by `sys.isotopes={'E','E','1H'}`.
- Lines 22-23: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 26-27: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 30-31: Sequence parameters; implemented by `parameters.fields=50e-6`.
- Lines 39-40: Simulation; implemented by `[yield,grid]=powder(spin_system,@rydmr_exp,parameters,'labframe')`.
- Lines 42-43: Preprocessing; implemented by `yield=cell2mat(yield)`.
- Lines 46-47: Plotting; implemented by `hull=get_hull(grid.betas,grid.gammas)`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=1`.
- Lines 16: computes `sys.isotopes` using `sys.isotopes={'E','E','1H'}`.
- Lines 17: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={2.0023 2.0023 0}`.
- Lines 18: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(3)`.
- Lines 19: computes `inter.coupling.matrix{1,3}` using `inter.coupling.matrix{1,3}=gauss2mhz([5e6 0 0`.
- Lines 23: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 24: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 27: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 31: computes `parameters.fields` using `parameters.fields=50e-6`.
- Lines 32: computes `parameters.rates` using `parameters.rates=2e6`.
- Lines 33: computes `parameters.electrons` using `parameters.electrons=[1 2]`.
- Lines 34: computes `parameters.grid` using `parameters.grid='leb_2ang_rank_71'`.
- Lines 35: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 36: computes `parameters.needs` using `parameters.needs={'zeeman_op'}`.
- Lines 37: computes `parameters.sum_up` using `parameters.sum_up=0`.
- Lines 40: computes `[yield,grid]` using `[yield,grid]=powder(spin_system,@rydmr_exp,parameters,'labframe')`.
- Lines 43: computes `yield` using `yield=cell2mat(yield)`.
- Lines 47: computes `hull` using `hull=get_hull(grid.betas,grid.gammas)`.

## Implementation structure

- Singlet yield anisotropy calculation for a radical pair
- using exponential recombination kinetics model.
- Calculation time: seconds
- Unit magnet (field sweep)
- System specification
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Preprocessing
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gauss2mhz()`, `create()`, `basis()`, `powder()`, `cell2mat()`, `get_hull()`, `kfigure()`, `trisurf()`.
