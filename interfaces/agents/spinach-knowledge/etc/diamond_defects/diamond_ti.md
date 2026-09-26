# etc/diamond_defects/diamond_ti.m

- Signature: `[sys,inter]=diamond_ti(parameters)`

## Purpose

Titanium-related defect spin system for diamond. Syntax: [sys,inter]=diamond_ti(parameters) Magnetic parameters from Nadolinny et al., Crystals 7, 237 (2017),

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- parameters is a structure with the following fields:
- .centre -'n3' or 'ok1'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .titanium -titanium isotope label, or 'none'
- .n_13c -number of reported 13C hyperfine couplings
- to include, from 0 to 2; applies only to OK1

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- Titanium-related defect spin system for diamond. Syntax:
- [sys,inter]=diamond_ti(parameters)
- Magnetic parameters from Nadolinny et al., Crystals 7, 237 (2017),
- parameters is a structure with the following fields:
- .centre -'n3' or 'ok1'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .titanium -titanium isotope label, or 'none'
- .n_13c -number of reported 13C hyperfine couplings
- to include, from 0 to 2; applies only to OK1
- sys -Spinach system specification structure
- inter -Spinach interaction specification structure
