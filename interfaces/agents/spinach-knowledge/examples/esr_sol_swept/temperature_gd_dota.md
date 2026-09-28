# examples/esr_sol_swept/temperature_gd_dota.m

- Signature: `temperature_gd_dota()`

## Purpose

Powder averaged W-band field-swept ESR spectrum of Gd(III) DOTA complex. Exact diagonalisation is used and a temperature dependence plot is produced. Calculation time: seconds.

## Physical / mathematical content

- The model is Gd(III) with isotropic g = 1.9918 and a zero-field-splitting tensor; it computes powder-averaged W-band field-swept ESR spectra at four temperatures.

## Numerical / algorithmic content

- For each temperature, the script rebuilds the spin system with exact diagonalisation in the Zeeman Hilbert-space formalism, then evaluates the field sweep on the same powder grid and plots the four spectra.

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
