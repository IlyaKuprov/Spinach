# etc/textbook/rlx_hfc.m

- Signature: `[r1,r2,rx]=rlx_hfc(B0,HFC,spins,tau_c)`

## Purpose

Redfield theory expressions for hyperfine relaxation and cross- relaxation rates, isotropic tumbling in liquid phase. Syntax: [r1,r2,rx]=rlx_hfc(B0,A,spins,tau_c)

## Physical / mathematical content

- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

## Parameters / inputs

- B0 -magnet field, Tesla
- A -3x3 hyperfine coupling tensor,
- not necessarily symmetric, rad/s
- spins -the spins involved, e.g. {'E','15N'},
- one of those must be an electron
- tau_c -rotational correlation time, seconds

## Outputs

- r1 -two longitudinal relaxation rates, Hz
- r2 -two transverse relaxation rates, Hz
- rx -longitudinal cross-relaxation rate, Hz

## Implementation structure

- Redfield theory expressions for hyperfine relaxation and cross-
- relaxation rates, isotropic tumbling in liquid phase. Syntax:
- [r1,r2,rx]=rlx_hfc(B0,A,spins,tau_c)
- B0 -magnet field, Tesla
- A -3x3 hyperfine coupling tensor,
- not necessarily symmetric, rad/s
- spins -the spins involved, e.g. {'E','15N'},
- one of those must be an electron
- tau_c -rotational correlation time, seconds
- r1 -two longitudinal relaxation rates, Hz
- r2 -two transverse relaxation rates, Hz
- rx -longitudinal cross-relaxation rate, Hz
