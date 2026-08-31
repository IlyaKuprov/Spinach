# examples/relaxation_theory/sle_esr_nitroxide_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/sle_esr_nitroxide_2.m`
- Signature: `sle_esr_nitroxide_2()`
- Total lines: 60

## Purpose

Slow motion regime simulation of an ESR spectrum of a nitroxide radical. Set to reproduce Figure 2 from the paper by Concilio et al.: Calculation time: seconds

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Magnet field; implemented by `sys.magnet=0.3343`.
- Lines 16-17: Isotopes; implemented by `sys.isotopes={'14N','E'}`.
- Lines 19-20: Coupling Matrices; implemented by `inter.coupling.matrix=cell(2)`.
- Lines 25-26: Zeeman Interactions; implemented by `inter.zeeman.matrix=cell(1, 2)`.
- Lines 31-32: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 35-36: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 39-40: SLE parameters; implemented by `parameters.max_rank=7`.
- Lines 53-54: SLE simulation; implemented by `spectrum_sle=gridfree(spin_system,@slowpass,parameters,'esr')`.
- Lines 56-57: SLE plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum_sle),parameters)`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=0.3343`.
- Lines 17: computes `sys.isotopes` using `sys.isotopes={'14N','E'}`.
- Lines 20: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(2)`.
- Lines 21: computes `inter.coupling.matrix{1,2}` using `inter.coupling.matrix{1,2}=1e6*gauss2mhz([1.000e+001 3.544e+000 1.170e+001`.
- Lines 26: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1, 2)`.
- Lines 27: computes `inter.zeeman.matrix{2}` using `inter.zeeman.matrix{2}=[2.0065794 -0.0007548 -0.0032848`.
- Lines 32: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 33: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 36: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 40: computes `parameters.max_rank` using `parameters.max_rank=7`.
- Lines 41: computes `parameters.tau_c` using `parameters.tau_c=17e-9`.
- Lines 42: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','E')`.
- Lines 43: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','E')`.
- Lines 44: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 45: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 46: computes `parameters.sweep` using `parameters.sweep=[-2.2e8 2e8]`.
- Lines 47: computes `parameters.npoints` using `parameters.npoints=1650`.
- Lines 48: computes `parameters.zerofill` using `parameters.zerofill=1650`.

## Implementation structure

- Slow motion regime simulation of an ESR spectrum of a nitroxide radical.
- Set to reproduce Figure 2 from the paper by Concilio et al.:
- Calculation time: seconds
- Magnet field
- Isotopes
- Coupling Matrices
- Zeeman Interactions
- Basis set
- Spinach housekeeping
- SLE parameters
- SLE simulation
- SLE plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gauss2mhz()`, `create()`, `basis()`, `state()`, `gridfree()`, `kfigure()`, `plot_1d()`.
