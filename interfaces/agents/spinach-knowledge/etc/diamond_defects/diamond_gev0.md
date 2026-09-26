# etc/diamond_defects/diamond_gev0.m

- Signature: `[sys,inter]=diamond_gev0(parameters)`

## Purpose

GeV0 spin system for diamond. Syntax: [sys,inter]=diamond_gev0(parameters) Magnetic parameters from Nadolinny et al., Phys. Status Solidi A 213, 2623 (2016), https://doi.org/10.1002/pssa.201600211

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- parameters is a structure with the following fields:
- .germanium -'73Ge', 'none', or another germanium isotope
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- GeV0 spin system for diamond. Syntax:
- [sys,inter]=diamond_gev0(parameters)
- Magnetic parameters from Nadolinny et al., Phys. Status Solidi A
- 213, 2623 (2016), https://doi.org/10.1002/pssa.201600211
- parameters is a structure with the following fields:
- .germanium -'73Ge', 'none', or another germanium isotope
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- sys -Spinach system specification structure
- inter -Spinach interaction specification structure
- Check input count
- Check consistency
