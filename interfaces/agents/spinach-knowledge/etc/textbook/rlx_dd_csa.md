# etc/textbook/rlx_dd_csa.m

## Signature

`[A,B,X]=rlx_dd_csa(B0,tau_c,isotopes,deltas,coords)`

## Purpose

Calculates Redfield relaxation and cross-relaxation rates for two spin-1/2 particles coupled by dipole–dipole (DD) and chemical-shift-anisotropy (CSA) interactions. The results include the CSA and DD contributions to each spin's longitudinal and transverse rates, the corresponding TROSY-component rates, and the longitudinal cross-relaxation rate between the spins.

## Model and calculation

The routine evaluates the CSA, dipolar, and CSA–DD cross-correlation contributions using the rotational correlation time and the interaction geometry. It reports broad and narrow TROSY-component transverse rates by adding or subtracting the absolute cross-correlation contribution from the CSA and DD contributions. The rotational correlation time is used in the spectral-density terms.

## Inputs

- `B0`: magnetic field in tesla.
- `tau_c`: rotational correlation time in seconds; it must be positive.
- `isotopes`: two isotope labels, for example `{'13C','19F'}`.
- `deltas`: cell array containing the two real, symmetric 3-by-3 chemical-shift tensors in ppm.
- `coords`: cell array containing the two 1-by-3 Cartesian coordinate vectors in angstroms.

## Outputs

- `A` and `B`: structures containing, for each spin, the CSA (`csa`), dipolar (`dd`), and total (`total`) contributions to longitudinal (`r1`) and transverse (`r2`) relaxation rates.
- `A.trosy` and `B.trosy`: `dd`, `csa`, and CSA–DD cross-correlation (`xc`) contributions, plus broad- and narrow-component transverse rates (`total_bro` and `total_nar`).
- `X`: longitudinal cross-relaxation rate between spins A and B.
