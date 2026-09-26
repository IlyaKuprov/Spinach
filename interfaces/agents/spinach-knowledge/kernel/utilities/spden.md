# kernel/utilities/spden.m

- Signature: `J=spden(L,D,omega)`

## Purpose

Lorentzian spectral density function for rotational diffusion at the user-specified frequency. Syntax: J=spden(L,D,omega)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

## Parameters / inputs

- L -spherical rank, use 2 for common NMR
- mechanisms such as dipolar relaxation
- D -rotational diffusion coefficient, s^{-1}
- omega -frequency, rad/s

## Outputs

- J -spectral density function value

## Implementation structure

- Lorentzian spectral density function for rotational
- diffusion at the user-specified frequency. Syntax:
- J=spden(L,D,omega)
- L -spherical rank, use 2 for common NMR
- mechanisms such as dipolar relaxation
- D -rotational diffusion coefficient, s^{-1}
- omega -frequency, rad/s
- J -spectral density function value
- Check consistency
- Get the correlation time
- Get the spectral density
- Consistency enforcement
