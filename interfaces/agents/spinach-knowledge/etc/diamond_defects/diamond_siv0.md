# etc/diamond_defects/diamond_siv0.m

- Signature: `[sys,inter]=diamond_siv0(parameters)`

## Purpose

SiV0 spin system for diamond. Syntax: [sys,inter]=diamond_siv0(parameters) Magnetic parameters from Edmonds et al., Phys. Rev. B 77, 245205 (2008), https://doi.org/10.1103/PhysRevB.77.245205

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- parameters is a structure with the following required fields:
- .silicon -'29Si', 'none', or another silicon isotope
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .n_13c -number of reported nearest-neighbour 13C
- hyperfine couplings, between 0 and 6

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- SiV0 spin system for diamond. Syntax:
- [sys,inter]=diamond_siv0(parameters)
- Magnetic parameters from Edmonds et al., Phys. Rev. B 77,
- 245205 (2008), https://doi.org/10.1103/PhysRevB.77.245205
- parameters is a structure with the following required fields:
- .silicon -'29Si', 'none', or another silicon isotope
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .n_13c -number of reported nearest-neighbour 13C
- hyperfine couplings, between 0 and 6
- sys -Spinach system specification structure
- inter -Spinach interaction specification structure
