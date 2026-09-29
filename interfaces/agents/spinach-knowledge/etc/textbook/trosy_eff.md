# etc/textbook/trosy_eff.m

## Use

eff=trosy_eff(B0,isotopes,xyz,csa) estimates how much of the first spin's chemical-shift-anisotropy (CSA) contribution to transverse linewidth is cancelled by dipole–dipole/CSA (DD–CSA) cross-correlation in a two-spin system.

## Inputs

- B0: finite, non-zero real magnetic field in tesla; the check does not require it to be positive.
- isotopes: cell array of two character isotope labels, both resolving to spin-1/2 isotopes (e.g. {'19F','13C'}).
- xyz: cell array of two real numeric three-element coordinate vectors with matching shape, giving the nuclei positions in Å.
- csa: real numeric 3-by-3 shielding or chemical-shift tensor for the first spin, in ppm. Its isotropic part is removed.

The nuclei must be at least 0.25 Å apart; closer coordinates are rejected. A non-zero second-rank anisotropic component of csa is required.

## Calculation and output

The routine obtains the dipolar tensor from xyz2dd, forms the field-scaled anisotropic Zeeman tensor from the traceless part of csa, and evaluates abs(X_DD_Z/DsqZ) using blprod and blinv. eff is a non-negative relative measure: the source describes eff = 1 as the case where the linewidth is purely dipolar. The implementation returns the absolute ratio directly and does not clamp it to an interval.

## Scope and source

Only two spin-1/2 isotopes are supported, and csa belongs to the first spin. Source: [implementation](https://github.com/IlyaKuprov/Spinach/blob/main/etc/textbook/trosy_eff.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=trosy_eff.m).
