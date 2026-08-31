# examples/spin_chemistry/singlet_yield_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/spin_chemistry/singlet_yield_2.m`
- Signature: `singlet_yield_2()`
- Total lines: 55

## Purpose

Liquid state magnetic field effect simulation on a radical pair with six equivalent nuclei using exponential recombi- nation kinetics model. Full S6 symmatry is used. Calculation time: seconds

## Physical / mathematical content

- Spin-chemistry examples. These scripts treat radical pairs, recombination channels, chemically induced dynamic nuclear polarisation, and magnetic-field effects. The theory combines spin-selective kinetics with singlet-triplet interconversion.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Unit magnet (field sweep); implemented by `sys.magnet=1`.
- Lines 16-17: System specification; implemented by `sys.isotopes ={'E','E','1H', '1H', '1H', '1H', '1H', '1H'}`.
- Lines 28-29: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 35-36: Fields and kinetics parameters; implemented by `parameters.rates=[0.176 0.880 1.76 3.52 8.8 17.6 35.2 52.8]*1e6`.
- Lines 42-43: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 46-47: Simulation; implemented by `M=liquid(spin_system,@rydmr_exp,parameters,'labframe')`.
- Lines 49-50: Plot the answer; implemented by `kfigure(); plot(parameters.fields,M,'r-')`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=1`.
- Lines 17: computes `sys.isotopes` using `sys.isotopes ={'E','E','1H', '1H', '1H', '1H', '1H', '1H'}`.
- Lines 18: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={2.002 2.002 0 0 0 0 0 0}`.
- Lines 19: computes `inter.coupling.scalar` using `inter.coupling.scalar=num2cell(mt2hz([0 0 0.295 0.295 0.295 0.295 0.295 0.295`.
- Lines 29: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 30: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 31: computes `bas.projections` using `bas.projections=0`.
- Lines 32: computes `bas.sym_spins` using `bas.sym_spins={[3 4 5 6 7 8]}`.
- Lines 33: computes `bas.sym_group` using `bas.sym_group={'S6'}`.
- Lines 36: computes `parameters.rates` using `parameters.rates=[0.176 0.880 1.76 3.52 8.8 17.6 35.2 52.8]*1e6`.
- Lines 37: computes `parameters.fields` using `parameters.fields=1e-3*(0:0.01:5)`.
- Lines 38: computes `parameters.electrons` using `parameters.electrons=[1 2]`.
- Lines 39: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 40: computes `parameters.needs` using `parameters.needs={'zeeman_op'}`.
- Lines 43: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 47: computes `M` using `M=liquid(spin_system,@rydmr_exp,parameters,'labframe')`.

## Implementation structure

- Liquid state magnetic field effect simulation on a radical
- pair with six equivalent nuclei using exponential recombi-
- nation kinetics model. Full S6 symmatry is used.
- Calculation time: seconds
- Unit magnet (field sweep)
- System specification
- Basis set
- Fields and kinetics parameters
- Spinach housekeeping
- Simulation
- Plot the answer

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `mt2hz()`, `create()`, `basis()`, `liquid()`, `kfigure()`, `kxlabel()`, `kylabel()`.
