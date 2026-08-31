# examples/dnp_mas/cross_effect_mas_enlev.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_mas/cross_effect_mas_enlev.m`
- Signature: `cross_effect_mas_enlev()`
- Total lines: 72

## Purpose

A MAS DNP simulation performed as described in Fred Mentink- Vigier's paper (Spinach rotation conventions are different): Energy level diagram as a function of the rotor phase. Calculation time: milliseconds

## Physical / mathematical content

- MAS DNP examples. These files model microwave-driven electron-nuclear polarisation transfer under magic-angle spinning, combining rotor-synchronised anisotropic interactions, relaxation, microwave irradiation, and powder/rotor averaging.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Magnet field; implemented by `sys.magnet=9.394`.
- Lines 17-18: Spin specification; implemented by `sys.isotopes={'E','E','1H'}`.
- Lines 20-21: Interactions; implemented by `inter.zeeman.eigs{1}=[2.0094 2.0060 2.0017]`.
- Lines 34-35: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 38-39: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 42-43: Stack generation parameters; implemented by `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 52-53: Stack generation; implemented by `H=rotor_stack(spin_system,parameters,'labframe')`.
- Lines 55-56: Stack diagonalization; implemented by `energies=zeros(size(H{1},1),numel(H))`.
- Lines 61-62: Plotting; implemented by `kfigure(); scale_figure([0.75 2.0])`.

### Control flow inferred from the code

- Line 57: `parfor` loop over `n=1:numel(H)`.

### Key state/data transformations

- Lines 15: computes `sys.magnet` using `sys.magnet=9.394`.
- Lines 18: computes `sys.isotopes` using `sys.isotopes={'E','E','1H'}`.
- Lines 21: computes `inter.zeeman.eigs{1}` using `inter.zeeman.eigs{1}=[2.0094 2.0060 2.0017]`.
- Lines 22: computes `inter.zeeman.euler{1}` using `inter.zeeman.euler{1}=[0.00 0.00 0.00]`.
- Lines 23: computes `inter.zeeman.eigs{2}` using `inter.zeeman.eigs{2}=[2.0094 2.0060 2.0017]`.
- Lines 24: computes `inter.zeeman.euler{2}` using `inter.zeeman.euler{2}=pi*[107 108 124]/180`.
- Lines 25: computes `inter.zeeman.eigs{3}` using `inter.zeeman.eigs{3}=[0.00 0.00 0.00]`.
- Lines 26: computes `inter.zeeman.euler{3}` using `inter.zeeman.euler{3}=[0.00 0.00 0.00]`.
- Lines 27: computes `inter.coupling.eigs` using `inter.coupling.eigs=cell(3,3)`.
- Lines 28: computes `inter.coupling.euler` using `inter.coupling.euler=cell(3,3)`.
- Lines 29: computes `inter.coupling.eigs{1,2}` using `inter.coupling.eigs{1,2}=[23.0e6 -11.5e6 -11.5e6]`.
- Lines 30: computes `inter.coupling.euler{1,2}` using `inter.coupling.euler{1,2}=pi*[0.00 135 0.00]/180`.
- Lines 31: computes `inter.coupling.eigs{1,3}` using `inter.coupling.eigs{1,3}=[1.5e6 -0.75e6 -0.75e6]`.
- Lines 32: computes `inter.coupling.euler{1,3}` using `inter.coupling.euler{1,3}=[0.00 0.00 0.00]`.
- Lines 35: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 36: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 39: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 43: computes `parameters.axis` using `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.

## Implementation structure

- A MAS DNP simulation performed as described in Fred Mentink-
- Vigier's paper (Spinach rotation conventions are different):
- Energy level diagram as a function of the rotor phase.
- Calculation time: milliseconds
- Magnet field
- Spin specification
- Interactions
- Basis set
- Spinach housekeeping
- Stack generation parameters
- Stack generation
- Stack diagonalization

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `rotor_stack()`, `energies()`, `kfigure()`, `scale_figure()`, `subplot()`, `kxlabel()`, `kylabel()`.
