# etc/diamond_defects/diamond_nv0_es.m

- Signature: `[sys,inter]=diamond_nv0_es(parameters)`

## Purpose

NV0 excited-state spin system for diamond. Syntax: [sys,inter]=diamond_nv0_es(parameters) Magnetic parameters from Felton et al., Phys. Rev. B 77, 081201 (2008), https://doi.org/10.1103/PhysRevB.77.081201

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- parameters is a structure with the following required fields:
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .nitrogen -'14N' or '15N'; 14N hyperfine couplings
- are scaled from 15N, with no NQI included
- because none is reported for this state

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- NV0 excited-state spin system for diamond. Syntax:
- [sys,inter]=diamond_nv0_es(parameters)
- Magnetic parameters from Felton et al., Phys. Rev. B 77,
- 081201 (2008), https://doi.org/10.1103/PhysRevB.77.081201
- parameters is a structure with the following required fields:
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .nitrogen -'14N' or '15N'; 14N hyperfine couplings
- are scaled from 15N, with no NQI included
- because none is reported for this state
- sys -Spinach system specification structure
- inter -Spinach interaction specification structure
