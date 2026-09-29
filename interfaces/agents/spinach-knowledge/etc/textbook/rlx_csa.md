# etc/textbook/rlx_csa.m

- MATLAB implementation: [etc/textbook/rlx_csa.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/textbook/rlx_csa.m)

**Signature:** `[r1,r2]=rlx_csa(B0,isotope,Z,tau_c)`

## Purpose

Calculates longitudinal and transverse Redfield relaxation rates for chemical-shift anisotropy (CSA), including the contributions associated with the antisymmetric part. The source notes that these rates do not depend on spin quantum number.

## Inputs

- `B0` — real scalar magnetic field in tesla.
- `isotope` — character-array isotope label, for example `'15N'`.
- `Z` — real 3-by-3 chemical-shift tensor in ppm.
- `tau_c` — positive real scalar second-rank rotational correlation time, documented as `1/(6D)`, in seconds.

The code checks that `B0` is a real numeric scalar, `tau_c` is positive, `isotope` is a character array, and `Z` is a real 3-by-3 matrix; it does not impose symmetry on `Z`.

## Calculation and outputs

The carrier-frequency expression is `omega=B0*(1+trace(1e-6*Z)/3)*spin(isotope)`. The tensor is converted from ppm by `1e-6`, then `blinv` supplies the invariants `Lsq` and `Dsq`. The implemented rates are

- `r1 = (1/2)Lsq*omega^2*tau_c/(1+9*tau_c^2*omega^2) + (2/15)Dsq*omega^2*tau_c/(1+tau_c^2*omega^2)`
- `r2 = (1/4)Lsq*omega^2*tau_c/(1+9*tau_c^2*omega^2) + (1/45)Dsq*omega^2*tau_c*(4+3/(1+tau_c^2*omega^2))`

Outputs `r1` and `r2` are respectively longitudinal and transverse rates in Hz.

## Reference

[Spinach Wiki: rlx_csa.m](https://spindynamics.org/wiki/index.php?title=rlx_csa.m)
