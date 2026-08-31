# etc/textbook/levelpop.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/textbook/levelpop.m`
- Signature: `[E,P,dP]=levelpop(isotope,field,temperature)`
- Total lines: 83

## Purpose

Equilibrium populations of the energy levels of a user-specified spin at the user-specified temperature. Energies are reported as fractions of kT at the temperature specified. Syntax: [E,P,dP]=levelpop(isotope,field,temperature)

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 36-37: Check consistency; implemented by `grumble(isotope,field,temperature)`.
- Lines 39-40: Fundamental constants; implemented by `h_bar=6.62607015e-34/(2*pi)`.
- Lines 43-44: Get the spin data; implemented by `[mg_ratio,multipl]=spin(isotope)`.
- Lines 46-47: Get Pauli matrices; implemented by `S=pauli(multipl)`.
- Lines 49-50: Build Zeeman Hamiltonian; implemented by `H=-mg_ratio*field*S.z; H=full(H)`.
- Lines 52-53: Get the energies; implemented by `E=h_bar*diag(H)/(k_bol*temperature)`.
- Lines 55-56: Get the populations; implemented by `P=exp(-h_bar*diag(H)/(k_bol*temperature))`.
- Lines 58-59: Normalise populations; implemented by `P=P/sum(P)`.
- Lines 61-62: Population differences; implemented by `dP=-diff(P)`.

### Key state/data transformations

- Lines 40: computes `h_bar` using `h_bar=6.62607015e-34/(2*pi)`.
- Lines 41: computes `k_bol` using `k_bol=1.380649e-23`.
- Lines 44: computes `[mg_ratio,multipl]` using `[mg_ratio,multipl]=spin(isotope)`.
- Lines 47: computes `S` using `S=pauli(multipl)`.
- Lines 50: computes `H` using `H=-mg_ratio*field*S.z; H=full(H)`.
- Lines 53: computes `E` using `E=h_bar*diag(H)/(k_bol*temperature)`.
- Lines 56: computes `P` using `P=exp(-h_bar*diag(H)/(k_bol*temperature))`.
- Lines 62: computes `dP` using `dP=-diff(P)`.

### Local helper functions

- Line 67: `grumble()` — `function grumble(isotope,field,temperature)`.
  - Representative operation: `if ~ischar(isotope)`.
  - Representative operation: `error('isotope must be a character string.')`.

## Parameters / inputs

- isotope -character string specifying the isotope.
- e.g. '1H', '13C', 'E', etc.
- field -primary magnet field in Tesla
- temperature -spin temperature, Kelvin

## Outputs

- E -vector of level energies, frac-
- tions of kT at the temperature
- specified
- P -vector of level populations
- dP -vector of population differences
- for adjacent levels
- Notes: the function is sensitive to the sign of the magnetogyric
- ratio -negative for electrons, positive for protons, etc.

## Implementation structure

- Equilibrium populations of the energy levels of a user-specified spin at
- the user-specified temperature. Energies are reported as fractions of kT
- at the temperature specified. Syntax:
- [E,P,dP]=levelpop(isotope,field,temperature)
- isotope -character string specifying the isotope.
- e.g. '1H', '13C', 'E', etc.
- field -primary magnet field in Tesla
- temperature -spin temperature, Kelvin
- E -vector of level energies, frac-
- tions of kT at the temperature
- specified
- P -vector of level populations

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `pauli()`, `diff()`, `ischar()`, `isscalar()`.
