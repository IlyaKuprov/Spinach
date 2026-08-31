# examples/spin_chemistry/singlet_yield_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/spin_chemistry/singlet_yield_1.m`
- Signature: `singlet_yield_1()`
- Total lines: 56

## Purpose

Liquid state magnetic field effect simulation on a radical pair with four nuclei using exponential recombination ki- netics model. Calculation time: seconds

## Physical / mathematical content

- Spin-chemistry examples. These scripts treat radical pairs, recombination channels, chemically induced dynamic nuclear polarisation, and magnetic-field effects. The theory combines spin-selective kinetics with singlet-triplet interconversion.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Unit magnet (field sweep); implemented by `sys.magnet=1`.
- Lines 16-17: System specification; implemented by `sys.isotopes ={'E','E','1H','1H','1H','1H'}`.
- Lines 26-27: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 33-34: Fields and kinetics parameters; implemented by `parameters.rates=[0.176 0.880 1.76 3.52 8.8 17.6 35.2 52.8]*1e6`.
- Lines 40-41: Disable ZTE; implemented by `sys.disable={'zte'}`.
- Lines 43-44: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 47-48: Simulation; implemented by `M=liquid(spin_system,@rydmr_exp,parameters,'labframe')`.
- Lines 50-51: Plot the answer; implemented by `kfigure(); plot(parameters.fields,M,'r-')`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=1`.
- Lines 17: computes `sys.isotopes` using `sys.isotopes ={'E','E','1H','1H','1H','1H'}`.
- Lines 18: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={2.002 2.002 0 0 0 0}`.
- Lines 19: computes `inter.coupling.scalar` using `inter.coupling.scalar = num2cell(mt2hz([0 0 0.195 0.195 0 0`.
- Lines 27: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 28: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 29: computes `bas.projections` using `bas.projections=0`.
- Lines 30: computes `bas.sym_spins` using `bas.sym_spins={[3 4]}`.
- Lines 31: computes `bas.sym_group` using `bas.sym_group={'S2'}`.
- Lines 34: computes `parameters.rates` using `parameters.rates=[0.176 0.880 1.76 3.52 8.8 17.6 35.2 52.8]*1e6`.
- Lines 35: computes `parameters.fields` using `parameters.fields=1e-3*(0:0.01:5)`.
- Lines 36: computes `parameters.electrons` using `parameters.electrons=[1 2]`.
- Lines 37: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 38: computes `parameters.needs` using `parameters.needs={'zeeman_op'}`.
- Lines 41: computes `sys.disable` using `sys.disable={'zte'}`.
- Lines 44: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 48: computes `M` using `M=liquid(spin_system,@rydmr_exp,parameters,'labframe')`.

## Implementation structure

- Liquid state magnetic field effect simulation on a radical
- pair with four nuclei using exponential recombination ki-
- netics model.
- Calculation time: seconds
- Unit magnet (field sweep)
- System specification
- Basis set
- Fields and kinetics parameters
- Disable ZTE
- Spinach housekeeping
- Simulation
- Plot the answer

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `mt2hz()`, `create()`, `basis()`, `liquid()`, `kfigure()`, `kxlabel()`, `kylabel()`.
