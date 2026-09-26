# etc/diamond_defects/diamond_vacancy.m

- Signature: `[sys,inter]=diamond_vacancy(parameters)`

## Purpose

Vacancy-family defect spin systems for diamond. Syntax: [sys,inter]=diamond_vacancy(parameters) R4/W6 parameters from Twitchen et al., Phys. Rev. B 59, 12900 (1999), https://doi.org/10.1103/PhysRevB.59.12900; W29 parameters from Kirui et al., Diam. Relat. Mater. 8, 1569 (1999), https://doi.org/10.1016/S0925-9635(99)00037-0; R5/O1/R6/R10/R11 parameters from Iakoubovskii and Stesmans, Phys. Rev. B 66, 045406 (2002), cr

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- parameters is a structure with the following fields:
- .centre -'r4_w6', 'w6', 'r4', 'w29', 'r5', 'o1',
- 'r6', 'r10', or 'r11'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- Vacancy-family defect spin systems for diamond. Syntax:
- [sys,inter]=diamond_vacancy(parameters)
- R4/W6 parameters from Twitchen et al., Phys. Rev. B 59, 12900
- (1999), https://doi.org/10.1103/PhysRevB.59.12900; W29
- parameters from Kirui et al., Diam. Relat. Mater. 8, 1569
- (1999), https://doi.org/10.1016/S0925-9635(99)00037-0;
- R5/O1/R6/R10/R11 parameters from Iakoubovskii and Stesmans,
- Phys. Rev. B 66, 045406 (2002),
- cross-checked against Ball, PhD thesis, OIST Graduate University
- (2021).
- parameters is a structure with the following fields:
- .centre -'r4_w6', 'w6', 'r4', 'w29', 'r5', 'o1',
