# etc/textbook/rlx_dip.m

## Signature

`[r1,r2,rx]=rlx_dip(B0,spins,dist,tau_c)`

## Purpose

Calculates Redfield longitudinal and transverse relaxation rates and longitudinal cross-relaxation for a pair of spins coupled by their through-space dipolar interaction in an isotropically tumbling liquid.

## Model and calculation

The routine constructs the dipolar interaction invariant from the inter-spin separation and combines it with each spin's spin-square factor. Spectral-density terms use the Zeeman frequencies, their sum and difference, and the zero-frequency contribution; the rotational diffusion coefficient is `1/(6*tau_c)`.

## Inputs

- `B0`: magnetic field in tesla.
- `spins`: two isotope labels, for example `{'1H','15N'}`.
- `dist`: inter-spin distance in angstroms.
- `tau_c`: positive rotational correlation time in seconds.

## Outputs

- `r1`: two longitudinal relaxation rates, in the order of the input spins, in Hz.
- `r2`: two transverse relaxation rates, in the order of the input spins, in Hz.
- `rx`: longitudinal cross-relaxation rate, in Hz.
