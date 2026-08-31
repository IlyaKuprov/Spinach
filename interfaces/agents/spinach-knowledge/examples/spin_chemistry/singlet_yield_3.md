# examples/spin_chemistry/singlet_yield_3.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/spin_chemistry/singlet_yield_3.m`
- Signature: `singlet_yield_3()`
- Total lines: 49

## Purpose

Figure 1 from the paper by Timmel, Till, Brocklehurst, McLauchlan and Hore: Calculation time: seconds

## Physical / mathematical content

- Spin-chemistry examples. These scripts treat radical pairs, recombination channels, chemically induced dynamic nuclear polarisation, and magnetic-field effects. The theory combines spin-selective kinetics with singlet-triplet interconversion.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Unit magnet (field sweep); implemented by `sys.magnet=1`.
- Lines 16-17: Spin system; implemented by `sys.isotopes={'E','E','1H'}`.
- Lines 19-20: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 23-24: Couplings; implemented by `inter.zeeman.scalar={2.0023 2.0023 0}`.
- Lines 28-29: Sequence parameters; implemented by `parameters.fields=linspace(0,3*20/1e4,200)`.
- Lines 36-37: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 40-41: Simulation; implemented by `M=liquid(spin_system,@rydmr_exp,parameters,'labframe')`.
- Lines 43-44: Plotting; implemented by `kfigure(); plot(linspace(0,3,200),M); kgrid`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=1`.
- Lines 17: computes `sys.isotopes` using `sys.isotopes={'E','E','1H'}`.
- Lines 20: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 21: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 24: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={2.0023 2.0023 0}`.
- Lines 25: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=gauss2mhz(20)*1e6`.
- Lines 26: computes `inter.coupling.scalar{3,3}` using `inter.coupling.scalar{3,3}=0`.
- Lines 29: computes `parameters.fields` using `parameters.fields=linspace(0,3*20/1e4,200)`.
- Lines 30-31: computes `parameters.rates` using `parameters.rates=2*pi*[0.005 0.02 0.05 0.1 0.15 0.2 0.3 0.5 2.0]*gauss2mhz(20)*1e6`.
- Lines 32: computes `parameters.electrons` using `parameters.electrons=[1 2]`.
- Lines 33: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 34: computes `parameters.needs` using `parameters.needs={'zeeman_op'}`.
- Lines 37: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 41: computes `M` using `M=liquid(spin_system,@rydmr_exp,parameters,'labframe')`.

## Implementation structure

- Figure 1 from the paper by Timmel, Till, Brocklehurst, McLauchlan
- and Hore:
- Calculation time: seconds
- Unit magnet (field sweep)
- Spin system
- Basis set
- Couplings
- Sequence parameters
- Spinach housekeeping
- Simulation
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gauss2mhz()`, `create()`, `basis()`, `liquid()`, `kfigure()`, `kylabel()`, `kxlabel()`.
