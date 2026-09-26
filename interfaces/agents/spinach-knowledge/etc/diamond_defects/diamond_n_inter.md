# etc/diamond_defects/diamond_n_inter.m

- Signature: `[sys,inter]=diamond_n_inter(parameters)`

## Purpose

Nitrogen interstitial spin system for diamond. Syntax: [sys,inter]=diamond_n_inter(parameters) Magnetic parameters from Felton et al., J. Phys. Condens. Matter 21, 364212 (2009), https://doi.org/10.1088/0953-8984/21/36/364212

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- parameters is a structure with the following fields:
- .centre -'war9' or 'war10'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .nitrogen -'14N' or '15N'

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- Nitrogen interstitial spin system for diamond. Syntax:
- [sys,inter]=diamond_n_inter(parameters)
- Magnetic parameters from Felton et al., J. Phys. Condens. Matter
- 21, 364212 (2009), https://doi.org/10.1088/0953-8984/21/36/364212
- parameters is a structure with the following fields:
- .centre -'war9' or 'war10'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .nitrogen -'14N' or '15N'
- sys -Spinach system specification structure
- inter -Spinach interaction specification structure
- Check input count
