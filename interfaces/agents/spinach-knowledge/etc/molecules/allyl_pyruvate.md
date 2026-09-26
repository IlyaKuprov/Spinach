# etc/molecules/allyl_pyruvate.m

- Signature: `[sys,inter]=allyl_pyruvate(spins)`

## Purpose

Spin system of allyl pyruvate. Isotropic chemical shifts and J-couplings determined by spectral fitting, coordinates and chemical shift anisotropies from DFT. Syntax: [sys,inter]=allyl_pyruvate(spins)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- spins -a cell array containing the isotopes
- to import, e.g. {'1H','13C'}

## Outputs

- sys, inter -Spinach data structures with the
- specification of the spin system
- Note: 13C-13C J-couplings are not provided -this spin system
- is for natural abundance 13C simulations only.

## Implementation structure

- Spin system of allyl pyruvate. Isotropic chemical shifts and
- J-couplings determined by spectral fitting, coordinates and
- chemical shift anisotropies from DFT. Syntax:
- [sys,inter]=allyl_pyruvate(spins)
- spins -a cell array containing the isotopes
- to import, e.g. {'1H','13C'}
- sys, inter -Spinach data structures with the
- specification of the spin system
- Note: 13C-13C J-couplings are not provided -this spin system
- is for natural abundance 13C simulations only.
- Check consistency
- Spin system
