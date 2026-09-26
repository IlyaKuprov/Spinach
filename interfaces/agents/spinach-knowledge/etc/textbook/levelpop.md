# etc/textbook/levelpop.m

- Signature: `[E,P,dP]=levelpop(isotope,field,temperature)`

## Purpose

Equilibrium populations of the energy levels of a user-specified spin at the user-specified temperature. Energies are reported as fractions of kT at the temperature specified. Syntax: [E,P,dP]=levelpop(isotope,field,temperature)

## Physical / mathematical content

## Numerical / algorithmic content

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
