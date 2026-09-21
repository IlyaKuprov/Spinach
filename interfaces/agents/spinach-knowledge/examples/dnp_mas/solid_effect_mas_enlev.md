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
