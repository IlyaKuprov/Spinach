# etc/diamond_defects/diamond_co.m

- Signature: `[sys,inter]=diamond_co(parameters)`

## Purpose

Cobalt-related defect spin system for diamond. Syntax: [sys,inter]=diamond_co(parameters) Magnetic parameters from Nadolinny et al., Crystals 7, 237 (2017),

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- parameters is a structure with the following fields:
- .centre -'o4' or 'nlo2'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- Cobalt-related defect spin system for diamond. Syntax:
- [sys,inter]=diamond_co(parameters)
- Magnetic parameters from Nadolinny et al., Crystals 7, 237 (2017),
- parameters is a structure with the following fields:
- .centre -'o4' or 'nlo2'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- sys -Spinach system specification structure
- inter -Spinach interaction specification structure
- Check input count
- Check consistency
- Set field-unit conversion constants
