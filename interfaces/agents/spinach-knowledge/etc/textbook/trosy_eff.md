# etc/textbook/trosy_eff.m

- Signature: `eff=trosy_eff(B0,isotopes,xyz,csa)`

## Purpose

TROSY efficiency in a two-spin system. Returns the extent of the cancellation of the CSA contribution to the trans- verse relaxation by the DD-CSA cross-correlation. Syntax: eff=trosy_eff(B0,isotopes,xyz,csa)

## Physical / mathematical content

- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

## Parameters / inputs

- B0 -magnet field, Tesla
- isotopes -a cell array with two character
- strings, e.g. {'19F','13C'} spe-
- cifying spin-1/2 isotopes
- xyz -a cell array with two Cartesian
- coordinate vectors in angstrom,
- giving the locations of the two
- nuclei
- csa -3x3 chemical shielding or chemi-
- cal shift (does not matter here)
- tensor of the first spin in ppm;
- its isotropic part will be drop-
- ped automatically

## Outputs

- eff -fraction of the CSA line width
- that is compensated by DD-CSA
- cross-correlation, 1 means that
- the line width is purely dipolar

## Implementation structure

- TROSY efficiency in a two-spin system. Returns the extent
- of the cancellation of the CSA contribution to the trans-
- verse relaxation by the DD-CSA cross-correlation. Syntax:
- eff=trosy_eff(B0,isotopes,xyz,csa)
- B0 -magnet field, Tesla
- isotopes -a cell array with two character
- strings, e.g. {'19F','13C'} spe-
- cifying spin-1/2 isotopes
- xyz -a cell array with two Cartesian
- coordinate vectors in angstrom,
- giving the locations of the two
- nuclei
