# etc/textbook/rlx_csa.m

## Signature

`[r1,r2]=rlx_csa(B0,isotope,Z,tau_c)`

## Purpose

Calculates longitudinal and transverse relaxation rates arising from chemical-shift anisotropy (CSA), using Redfield theory for isotropic rotational motion. The expression includes contributions associated with the antisymmetric part of the CSA tensor.

## Model and calculation

The routine obtains the carrier frequency from the field, isotope, and isotropic part of the chemical-shift tensor. It then evaluates the Blicharski invariants of the tensor and combines them with Lorentzian spectral-density factors at the relevant frequencies. The isotropic rotational correlation time is supplied as the second-rank correlation time, `tau_c=1/(6D)`.

## Inputs

- `B0`: magnetic field in tesla.
- `isotope`: spin isotope label, for example `'15N'`.
- `Z`: real 3-by-3 chemical-shift tensor in ppm.
- `tau_c`: second-rank rotational correlation time in seconds.

## Outputs

- `r1`: longitudinal relaxation rate in Hz.
- `r2`: transverse relaxation rate in Hz.

The expressions do not depend on the spin quantum number.
