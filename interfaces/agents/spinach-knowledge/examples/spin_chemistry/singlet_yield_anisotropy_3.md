# examples/spin_chemistry/singlet_yield_anisotropy_3.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/spin_chemistry/singlet_yield_anisotropy_3.m`
- Signature: `singlet_yield_anisotropy_3()`
- Total lines: 72

## Purpose

Singlet yield anisotropy calculation for a model radical pair reaction, Haberkorn recombination model. Run time: minutes on NVidia Titan V card, hours on CPU.

## Physical / mathematical content

- Spin-chemistry examples. These scripts treat radical pairs, recombination channels, chemically induced dynamic nuclear polarisation, and magnetic-field effects. The theory combines spin-selective kinetics with singlet-triplet interconversion.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Earth field; implemented by `sys.magnet=50e-6`.
- Lines 15-16: Isotopes; implemented by `sys.isotopes={'E','E','14N','14N','1H','1H','1H'}`.
- Lines 18-19: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 23-24: Hyperfine coupling tensors; implemented by `inter.coupling.matrix=cell(7)`.
- Lines 41-42: Zeeman interactions; implemented by `inter.zeeman.scalar={2.0023 2.0025 0 0 0 0 0}`.
- Lines 44-45: Kinetics parameters; implemented by `inter.chem.rp_theory='haberkorn'`.
- Lines 49-50: Sequence parameters; implemented by `parameters.grid='leb_1ang_rank_63'`.
- Lines 56-60: Enable GPU arithmetic sys.enable={'gpu'};; implemented by `spin_system=create(sys,inter)`.
- Lines 59-60: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 63-64: Run a simulation; implemented by `[yields,grid]=powder(spin_system,@rydmr,parameters,'labframe')`.
- Lines 66-67: Do the plotting; implemented by `kfigure(); plot(grid.betas,cell2mat(yields))`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=50e-6`.
- Lines 16: computes `sys.isotopes` using `sys.isotopes={'E','E','14N','14N','1H','1H','1H'}`.
- Lines 19: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 20: computes `bas.approximation` using `bas.approximation='IK-0'`.
- Lines 21: computes `bas.level` using `bas.level=5`.
- Lines 24: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(7)`.
- Lines 25: computes `inter.coupling.matrix{1,3}` using `inter.coupling.matrix{1,3}=mt2hz([-0.0989 0.0039 0.0000`.
- Lines 28: computes `inter.coupling.matrix{2,4}` using `inter.coupling.matrix{2,4}=mt2hz([-0.0336 0.0924 -0.1354`.
- Lines 31: computes `inter.coupling.matrix{2,5}` using `inter.coupling.matrix{2,5}=mt2hz([-0.9920 -0.2091 -0.2003`.
- Lines 34: computes `inter.coupling.matrix{2,6}` using `inter.coupling.matrix{2,6}=mt2hz([-0.9920 -0.2091 -0.2003`.
- Lines 37: computes `inter.coupling.matrix{2,7}` using `inter.coupling.matrix{2,7}=mt2hz([-0.9920 -0.2091 -0.2003`.
- Lines 42: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={2.0023 2.0025 0 0 0 0 0}`.
- Lines 45: computes `inter.chem.rp_theory` using `inter.chem.rp_theory='haberkorn'`.
- Lines 46: computes `inter.chem.rp_electrons` using `inter.chem.rp_electrons=[1 2]`.
- Lines 47: computes `inter.chem.rp_rates` using `inter.chem.rp_rates=[1e6 1e6]`.
- Lines 50: computes `parameters.grid` using `parameters.grid='leb_1ang_rank_63'`.
- Lines 51: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 52: computes `parameters.tol` using `parameters.tol=1e-2`.

## Implementation structure

- Singlet yield anisotropy calculation for a model radical pair
- reaction, Haberkorn recombination model.
- Run time: minutes on NVidia Titan V card, hours on CPU.
- Earth field
- Isotopes
- Basis set
- Hyperfine coupling tensors
- Zeeman interactions
- Kinetics parameters
- Sequence parameters
- Enable GPU arithmetic
- sys.enable={'gpu'};

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `mt2hz()`, `create()`, `basis()`, `powder()`, `kfigure()`, `cell2mat()`, `kxlabel()`, `kylabel()`.
