# etc/textbook/trosy_eff.m

## Signature

`eff=trosy_eff(B0,isotopes,xyz,csa)`

## Purpose

Estimates how much of the CSA contribution to transverse linewidth is cancelled by dipole–dipole/CSA (DD–CSA) cross-correlation for two spin-1/2 nuclei. The result is the magnitude of the cross-correlation invariant divided by the second-rank CSA invariant; the function describes `eff=1` as the limit in which the linewidth is purely dipolar.

## Model and calculation

The routine constructs the dipolar coupling tensor from the two coordinates and isotope labels. It removes the isotropic part of the first spin's chemical-shift tensor, scales its anisotropic part by the magnetic field and the spin gyromagnetic ratio, and evaluates the CSA invariant and its cross-invariant with the dipolar tensor.

## Inputs

- `B0`: magnetic field in tesla.
- `isotopes`: cell array of two spin-1/2 isotope labels, for example `{'19F','13C'}`.
- `xyz`: cell array containing the two three-element Cartesian nuclear coordinates in angstroms.
- `csa`: real 3-by-3 chemical-shift or shielding tensor for the first spin, in ppm. Its isotropic part is removed.

## Output

- `eff`: dimensionless magnitude ratio describing the CSA linewidth cancellation by DD–CSA cross-correlation.
