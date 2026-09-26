# etc/textbook/rlx_dip.m

- Signature: `[r1,r2,rx]=rlx_dip(B0,spins,dist,tau_c)`

## Purpose

Redfield theory expressions for dipolar relaxation and cross- relaxation rates, isotropic tumbling in liquid phase. Syntax: [r1,r2,rx]=rlx_dip(B0,spins,dist,tau_c)

## Physical / mathematical content

- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

## Parameters / inputs

- B0 -magnet field, Tesla
- spins -the spins involved, e.g. {'1H','15N'}
- dist -inter-spin distance, Angstrom
- tau_c -rotational correlation time, seconds

## Outputs

- r1 -two longitudinal relaxation rates, Hz
- r2 -two transverse relaxation rates, Hz
- rx -longitudinal cross-relaxation rate, Hz

## Implementation structure

- Redfield theory expressions for dipolar relaxation and cross-
- relaxation rates, isotropic tumbling in liquid phase. Syntax:
- [r1,r2,rx]=rlx_dip(B0,spins,dist,tau_c)
- B0 -magnet field, Tesla
- spins -the spins involved, e.g. {'1H','15N'}
- dist -inter-spin distance, Angstrom
- tau_c -rotational correlation time, seconds
- r1 -two longitudinal relaxation rates, Hz
- r2 -two transverse relaxation rates, Hz
- rx -longitudinal cross-relaxation rate, Hz
- Check consistency
- Blicharsky invariant and rotational diffusion coefficient
