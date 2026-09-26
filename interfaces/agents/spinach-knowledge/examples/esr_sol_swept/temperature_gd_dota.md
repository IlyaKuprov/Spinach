# examples/esr_sol_swept/temperature_gd_dota.m

- Signature: `temperature_gd_dota()`

## Purpose

Powder averaged W-band field-swept ESR spectrum of Gd(III) DOTA complex. Exact diagonalisation is used and a tempera- ture dependence plot is produced. Calculation time: seconds.

## Physical / mathematical content

- Field-swept ESR examples. These files emphasise resonance-field finding, powder averaging, anisotropic g and hyperfine tensors, and intensity accumulation over orientation manifolds.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Powder averaged W-band field-swept ESR spectrum of Gd(III)
- DOTA complex. Exact diagonalisation is used and a tempera-
- ture dependence plot is produced.
- Calculation time: seconds.
- Isotopes
- Magnet field (must be 1)
- Properties
- Basis set
- Get figure going
- Temperatures
- Loop over temperatures
- Set the temperature
