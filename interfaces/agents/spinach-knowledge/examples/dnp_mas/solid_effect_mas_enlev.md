# examples/dnp_mas/solid_effect_mas_enlev.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_mas/solid_effect_mas_enlev.m`
- Signature: `solid_effect_mas_enlev()`
- Total lines: 62

## Purpose

A MAS DNP simulation performed as described in Fred Mentink- Vigier's paper (Spinach rotation conventions are different): Energy level diagram as a function of the rotor phase. Calculation time: milliseconds

## Physical / mathematical content

- MAS DNP examples. These files model microwave-driven electron-nuclear polarisation transfer under magic-angle spinning, combining rotor-synchronised anisotropic interactions, relaxation, microwave irradiation, and powder/rotor averaging.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Magnet field; implemented by `sys.magnet=9.403`.
- Lines 17-18: Spin specification; implemented by `sys.isotopes={'E','1H'}`.
- Lines 20-21: Interactions; implemented by `inter.zeeman.eigs{1}=[2.00614 2.00194 2.00988]`.
- Lines 28-29: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 32-33: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 36-37: Stack generation parameters; implemented by `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 45-46: Stack generation; implemented by `H=rotor_stack(spin_system,parameters,'esr')`.
- Lines 48-49: Stack diagonalization; implemented by `energies=zeros(size(H{1},1),numel(H))`.
- Lines 54-55: Plotting; implemented by `phi_axis=linspace(0,2*pi,numel(H)+1)'`.

### Control flow inferred from the code

- Line 50: `parfor` loop over `n=1:numel(H)`.

### Key state/data transformations

- Lines 15: computes `sys.magnet` using `sys.magnet=9.403`.
- Lines 18: computes `sys.isotopes` using `sys.isotopes={'E','1H'}`.
- Lines 21: computes `inter.zeeman.eigs{1}` using `inter.zeeman.eigs{1}=[2.00614 2.00194 2.00988]`.
- Lines 22: computes `inter.zeeman.euler{1}` using `inter.zeeman.euler{1}=pi*[253.6 105.1 123.8]/180`.
- Lines 23: computes `inter.zeeman.eigs{2}` using `inter.zeeman.eigs{2}=[0.00 0.00 0.00]`.
- Lines 24: computes `inter.zeeman.euler{2}` using `inter.zeeman.euler{2}=[0.00 0.00 0.00]`.
- Lines 25-26: computes `inter.coordinates` using `inter.coordinates={[0.00 0.00 0.00], [0.00 0.00 3.00]}`.
- Lines 29: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 30: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 33: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 37: computes `parameters.axis` using `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 38: computes `parameters.rframes` using `parameters.rframes={}`.
- Lines 39: computes `parameters.orientation` using `parameters.orientation=[0 0 0]`.
- Lines 40: computes `parameters.spins` using `parameters.spins={'E','1H'}`.
- Lines 41: computes `parameters.masframe` using `parameters.masframe='magnet'`.
- Lines 42: computes `parameters.offset` using `parameters.offset=[0 0]`.
- Lines 43: computes `parameters.max_rank` using `parameters.max_rank=200`.
- Lines 46: computes `H` using `H=rotor_stack(spin_system,parameters,'esr')`.

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

- Called routines detected from the main body: `create()`, `basis()`, `rotor_stack()`, `energies()`, `phi_axis()`, `kfigure()`, `kxlabel()`, `kylabel()`.
