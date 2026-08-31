# examples/spin_chemistry/singlet_yield_4.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/spin_chemistry/singlet_yield_4.m`
- Signature: `singlet_yield_4()`
- Total lines: 54

## Purpose

Figure 3 from the paper by Till, Timmel, Brocklehurst and Hore: Note: the original paper only uses electron Zeeman operators for the field sweep, and therefore misses the effects associa- ted with the rise in the nuclear Zeeman interaction on the high field side of the resulting plot. Calculation time: seconds

## Physical / mathematical content

- Spin-chemistry examples. These scripts treat radical pairs, recombination channels, chemically induced dynamic nuclear polarisation, and magnetic-field effects. The theory combines spin-selective kinetics with singlet-triplet interconversion.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 18-19: Unit magnet (field sweep); implemented by `sys.magnet=1`.
- Lines 21-22: Spin system; implemented by `sys.isotopes={'E','E','1H','1H'}`.
- Lines 24-25: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 28-29: Couplings; implemented by `inter.zeeman.scalar={2.0023 2.0044 0 0}`.
- Lines 34-35: Sequence parameters; implemented by `parameters.fields=1e-3*10.^linspace(-5,3,2000)`.
- Lines 41-42: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 45-46: Simulation; implemented by `M=liquid(spin_system,@rydmr_exp,parameters,'labframe')`.
- Lines 48-49: Plotting; implemented by `kfigure(); plot(linspace(-5,3,2000),M); kgrid`.

### Key state/data transformations

- Lines 19: computes `sys.magnet` using `sys.magnet=1`.
- Lines 22: computes `sys.isotopes` using `sys.isotopes={'E','E','1H','1H'}`.
- Lines 25: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 26: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 29: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={2.0023 2.0044 0 0}`.
- Lines 30: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=gauss2mhz(35)*1e6`.
- Lines 31: computes `inter.coupling.scalar{1,4}` using `inter.coupling.scalar{1,4}=gauss2mhz(30)*1e6`.
- Lines 32: computes `inter.coupling.scalar{4,4}` using `inter.coupling.scalar{4,4}=0`.
- Lines 35: computes `parameters.fields` using `parameters.fields=1e-3*10.^linspace(-5,3,2000)`.
- Lines 36: computes `parameters.rates` using `parameters.rates=[0.1 1.0 10.0 100.0 1000.0]*1e6`.
- Lines 37: computes `parameters.electrons` using `parameters.electrons=[1 2]`.
- Lines 38: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 39: computes `parameters.needs` using `parameters.needs={'zeeman_op'}`.
- Lines 42: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 46: computes `M` using `M=liquid(spin_system,@rydmr_exp,parameters,'labframe')`.

## Implementation structure

- Figure 3 from the paper by Till, Timmel, Brocklehurst and Hore:
- Note: the original paper only uses electron Zeeman operators for
- the field sweep, and therefore misses the effects associa-
- ted with the rise in the nuclear Zeeman interaction on the
- high field side of the resulting plot.
- Calculation time: seconds
- Unit magnet (field sweep)
- Spin system
- Basis set
- Couplings
- Sequence parameters
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gauss2mhz()`, `create()`, `basis()`, `liquid()`, `kfigure()`, `kylabel()`, `kxlabel()`.
