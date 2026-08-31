# examples/spin_chemistry/singlet_yield_5.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/spin_chemistry/singlet_yield_5.m`
- Signature: `singlet_yield_5()`
- Total lines: 47

## Purpose

Figure 5 from the paper by Timmel, Till, Brocklehurst, McLauchlan and Hore: Calculation time: seconds

## Physical / mathematical content

- Spin-chemistry examples. These scripts treat radical pairs, recombination channels, chemically induced dynamic nuclear polarisation, and magnetic-field effects. The theory combines spin-selective kinetics with singlet-triplet interconversion.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Unit magnet (field sweep); implemented by `sys.magnet=1`.
- Lines 16-17: System specification; implemented by `sys.isotopes={'E','E','1H','1H'}`.
- Lines 23-24: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 27-28: Magnetic field and kinetics; implemented by `parameters.fields=linspace(0,0.5*20/1e4,200)`.
- Lines 34-35: Spinach run; implemented by `spin_system=create(sys,inter)`.
- Lines 38-39: Simulation; implemented by `M=liquid(spin_system,@rydmr_exp,parameters,'labframe')`.
- Lines 41-42: Plotting; implemented by `kfigure(); plot(linspace(0,0.5,200),M); kgrid`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=1`.
- Lines 17: computes `sys.isotopes` using `sys.isotopes={'E','E','1H','1H'}`.
- Lines 18: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={2.0023 2.0023 0 0}`.
- Lines 19: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=gauss2mhz(20)*1e6`.
- Lines 20: computes `inter.coupling.scalar{2,4}` using `inter.coupling.scalar{2,4}=gauss2mhz(20)*1e6`.
- Lines 21: computes `inter.coupling.scalar{4,4}` using `inter.coupling.scalar{4,4}=0`.
- Lines 24: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 25: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 28: computes `parameters.fields` using `parameters.fields=linspace(0,0.5*20/1e4,200)`.
- Lines 29: computes `parameters.rates` using `parameters.rates=2*pi*[1e-4 1e-3 0.01 0.02 0.05 0.1 0.2]*gauss2mhz(20)*1e6`.
- Lines 30: computes `parameters.electrons` using `parameters.electrons=[1 2]`.
- Lines 31: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 32: computes `parameters.needs` using `parameters.needs={'zeeman_op'}`.
- Lines 35: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 39: computes `M` using `M=liquid(spin_system,@rydmr_exp,parameters,'labframe')`.

## Implementation structure

- Figure 5 from the paper by Timmel, Till, Brocklehurst, McLauchlan
- and Hore:
- Calculation time: seconds
- Unit magnet (field sweep)
- System specification
- Basis set
- Magnetic field and kinetics
- Spinach run
- Simulation
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gauss2mhz()`, `create()`, `basis()`, `liquid()`, `kfigure()`, `kylabel()`, `kxlabel()`.
