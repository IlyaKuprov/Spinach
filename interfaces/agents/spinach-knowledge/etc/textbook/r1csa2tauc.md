# etc/textbook/r1csa2tauc.m

- MATLAB implementation: [etc/textbook/r1csa2tauc.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/textbook/r1csa2tauc.m)

- Signature: `tauc=r1csa2tauc(R1,del_sq,B0,isotope)`

## Purpose

Returns the two rotational correlation-time candidates compatible with a longitudinal CSA relaxation rate under the quadratic relation implemented here.

## Inputs

All four arguments are required; there are no defaults.

- R1 — positive real numeric scalar longitudinal relaxation rate, documented in Hz.
- del_sq — positive real numeric scalar second-rank invariant of the CSA. It is the dimensionless squared invariant of fractional shielding: pass the second output of [blinv.m](https://spindynamics.org/wiki/index.php?title=blinv.m) applied to `1e-6*Z_ppm`, not to the unconverted ppm tensor.
- B0 — real numeric scalar magnetic field, in tesla. The input check does not require it to be positive or nonzero, but the calculation requires a nonzero Zeeman frequency.
- isotope — character array naming the isotope, for example '1H'.

## Calculation and output

The code sets omega = -B0*spin(isotope), evaluates the plus-square-root quadratic candidate first in tauc(2) to reduce cancellation, then sets tauc(1) = 1/(omega^2*tauc(2)). The returned two-element row vector is therefore [tauc(1), tauc(2)]: the reciprocal-derived candidate is component 1 and the plus-branch candidate is component 2. This corrects the prior page's claim that the returned vector is ordered with the larger candidate first; the source computes the plus branch first but stores it at index 2. The source documents tauc in seconds.

The discriminant is del_sq^2*omega^4 - 225*omega^2*R1^2. The source raises its “no real solutions” error only when both computed components are non-real; it does not require both candidates individually to be real before returning.

## Source

[Spinach Wiki: r1csa2tauc.m](https://spindynamics.org/wiki/index.php?title=r1csa2tauc.m).
