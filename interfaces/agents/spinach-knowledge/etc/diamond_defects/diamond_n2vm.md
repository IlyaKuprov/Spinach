# etc/diamond_defects/diamond_n2vm.m

- Signature: `[sys,inter]=diamond_n2vm(parameters)`

## Purpose

N2V-spin system for diamond. Syntax: [sys,inter]=diamond_n2vm(parameters) Magnetic parameters from Green et al., Phys. Rev. B 92, 165204 (2015), https://doi.org/10.1103/PhysRevB.92.165204

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- parameters is a structure with the following required fields:
- .nitrogen -'14N' or '15N'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .include_13c -include reported 13C hyperfine couplings,
- true or false

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- N2V-spin system for diamond. Syntax:
- [sys,inter]=diamond_n2vm(parameters)
- Magnetic parameters from Green et al., Phys. Rev. B 92,
- 165204 (2015), https://doi.org/10.1103/PhysRevB.92.165204
- parameters is a structure with the following required fields:
- .nitrogen -'14N' or '15N'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .include_13c -include reported 13C hyperfine couplings,
- true or false
- sys -Spinach system specification structure
- inter -Spinach interaction specification structure
