# etc/diamond_defects/diamond_p.m

- Signature: `[sys,inter]=diamond_p(parameters)`

## Purpose

Phosphorus-related defect spin system for diamond. Syntax: [sys,inter]=diamond_p(parameters) Magnetic parameters from Nadolinny et al., Crystals 7, 237 (2017),

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- parameters is a structure with the following fields:
- .centre -'ma1', 'np1', 'np2', 'np3', 'np4', 'np5',
- 'np6', 'np8', or 'np9'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .include_13c -include the reported 13C hyperfine coupling;
- applies only to MA1 and defaults to false

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- Phosphorus-related defect spin system for diamond. Syntax:
- [sys,inter]=diamond_p(parameters)
- Magnetic parameters from Nadolinny et al., Crystals 7, 237 (2017),
- parameters is a structure with the following fields:
- .centre -'ma1', 'np1', 'np2', 'np3', 'np4', 'np5',
- 'np6', 'np8', or 'np9'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .include_13c -include the reported 13C hyperfine coupling;
- applies only to MA1 and defaults to false
- sys -Spinach system specification structure
- inter -Spinach interaction specification structure
