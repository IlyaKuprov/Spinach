# examples/relaxation_theory/sle_solid_limit.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/sle_solid_limit.m`
- Signature: `sle_solid_limit()`
- Total lines: 73

## Purpose

Solid limit of Stochastic Liouville equation formalism. Calculation time: hours

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Magnet field; implemented by `sys.magnet=0.3343`.
- Lines 12-13: Isotopes; implemented by `sys.isotopes={'14N','E'}`.
- Lines 15-16: Coupling Matrices; implemented by `inter.coupling.matrix=cell(2)`.
- Lines 21-22: Zeeman Interactions; implemented by `inter.zeeman.matrix=cell(1, 2)`.
- Lines 27-28: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 31-32: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 35-36: SLE parameters; implemented by `parameters.rho0=state(spin_system,'L+','E')`.
- Lines 47-48: Ranks and correlation times; implemented by `ranks=[3 7 15 30 ]`.
- Lines 51-52: Start a figure; implemented by `kfigure(); scale_figure([3.0 1.0])`.
- Lines 54-55: Loop over correlation times; implemented by `for n=1:numel(ranks)`.
- Lines 57-58: Set the parameters; implemented by `parameters.max_rank=ranks(n)`.
- Lines 61-62: SLE simulation; implemented by `spectrum_sle=gridfree(spin_system,@slowpass,parameters,'esr')`.
- Lines 64-65: SLE plotting; implemented by `subplot(1,numel(ranks),n)`.

### Control flow inferred from the code

- Line 55: `for` loop over `n=1:numel(ranks)`.

### Key state/data transformations

- Lines 10: computes `sys.magnet` using `sys.magnet=0.3343`.
- Lines 13: computes `sys.isotopes` using `sys.isotopes={'14N','E'}`.
- Lines 16: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(2)`.
- Lines 17: computes `inter.coupling.matrix{1,2}` using `inter.coupling.matrix{1,2}=1e6*gauss2mhz([1.000e+001 3.544e+000 1.170e+001`.
- Lines 22: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1, 2)`.
- Lines 23: computes `inter.zeeman.matrix{2}` using `inter.zeeman.matrix{2}=[2.0065794 -0.0007548 -0.0032848`.
- Lines 28: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 29: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 32: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 36: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','E')`.
- Lines 37: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','E')`.
- Lines 38: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 39: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 40: computes `parameters.sweep` using `parameters.sweep=[-2.2e8 2e8]`.
- Lines 41: computes `parameters.npoints` using `parameters.npoints=240`.
- Lines 42: computes `parameters.zerofill` using `parameters.zerofill=240`.
- Lines 43: computes `parameters.axis_units` using `parameters.axis_units='GHz-labframe'`.
- Lines 44: computes `parameters.invert_axis` using `parameters.invert_axis=1`.

## Implementation structure

- Solid limit of Stochastic Liouville equation formalism.
- Calculation time: hours
- Magnet field
- Isotopes
- Coupling Matrices
- Zeeman Interactions
- Basis set
- Spinach housekeeping
- SLE parameters
- Ranks and correlation times
- Start a figure
- Loop over correlation times

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gauss2mhz()`, `create()`, `basis()`, `state()`, `kfigure()`, `scale_figure()`, `ranks()`, `tau_c()`, `gridfree()`, `subplot()`, `plot_1d()`, `ktitle()`, `num2str()`, `log10()`, `kxlabel()`.
