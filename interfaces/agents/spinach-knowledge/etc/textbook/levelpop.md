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
