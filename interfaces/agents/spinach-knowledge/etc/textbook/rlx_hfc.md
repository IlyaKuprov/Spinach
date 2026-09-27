# etc/textbook/rlx_hfc.m

## Signature

`[r1,r2,rx]=rlx_hfc(B0,HFC,spins,tau_c)`

## Purpose

Calculates Redfield longitudinal and transverse relaxation rates and longitudinal cross-relaxation for a pair of spins coupled by a hyperfine interaction under isotropic tumbling in a liquid. One of the two spins must be an electron.

## Model and calculation

The routine obtains the rank-1 and rank-2 Blicharski invariants of the hyperfine tensor and combines both contributions with spectral densities at the spins' Zeeman frequencies, their sum and difference, and zero frequency. Spin-square factors account for the two spin quantum numbers; the rotational diffusion coefficient is `1/(6*tau_c)`.

## Inputs

- `B0`: magnetic field in tesla.
- `HFC`: real 3-by-3 hyperfine coupling tensor in radians per second; it need not be symmetric.
- `spins`: two isotope labels, one identifying an electron; for example, `{'E','15N'}`.
- `tau_c`: positive rotational correlation time in seconds.

## Outputs

- `r1`: two longitudinal relaxation rates, in the order of the input spins, in Hz.
- `r2`: two transverse relaxation rates, in the order of the input spins, in Hz.
- `rx`: longitudinal cross-relaxation rate, in Hz.
