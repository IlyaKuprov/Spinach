# etc/diamond_defects/diamond_p1.m

- Signature: `[sys,inter]=diamond_p1(parameters)`

## Purpose

P1 centre spin system for diamond. Syntax: [sys,inter]=diamond_p1(parameters) Magnetic parameters from: Nir-Arad et al. Phys. Chem. Chem. Phys. 26, 27633 (2024), <https://doi.org/10.1039/d4cp03055a>, and Smith et al. Phys. Rev. 115, 1546 (1959),

## Physical / mathematical content

- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

## Parameters / inputs

- a structure (parameters.*) with the following fields:
- .orientation -'111', '110', or '100' crystal
- plane normal aligned with the
- magnetic field, def. is '111'
- .nitrogen -'14N' or '15N', default is '14N'

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- P1 centre spin system for diamond. Syntax:
- [sys,inter]=diamond_p1(parameters)
- Magnetic parameters from: Nir-Arad et al. Phys. Chem. Chem. Phys. 26,
- 27633 (2024), <https://doi.org/10.1039/d4cp03055a>, and
- Smith et al. Phys. Rev. 115, 1546 (1959),
- a structure (parameters.*) with the following fields:
- .orientation -'111', '110', or '100' crystal
- plane normal aligned with the
- magnetic field, def. is '111'
- .nitrogen -'14N' or '15N', default is '14N'
- sys -Spinach system specification structure
- inter -Spinach interaction specification structure
