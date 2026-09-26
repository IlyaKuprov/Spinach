# etc/diamond_defects/diamond_nvm_gs.m

- Signature: `[sys,inter]=diamond_nvm_gs(parameters)`

## Purpose

NV centre ground state spin system for diamond. Syntax: [sys,inter]=diamond_nvm_gs(parameters) Magnetic parameters from: S. Felton et al., Phys. Rev. B 79, 075203 (2009),

## Physical / mathematical content

- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

## Parameters / inputs

- the following is needed in the parameters.* structure:
- .orientation -'111', '110', or '100' crystal
- plane normal aligned with the
- magnetic field, def. is '111'
- .nitrogen -'14N' or '15N', default is '14N'

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- NV centre ground state spin system for diamond. Syntax:
- [sys,inter]=diamond_nvm_gs(parameters)
- Magnetic parameters from:
- S. Felton et al., Phys. Rev. B 79, 075203 (2009),
- the following is needed in the parameters.* structure:
- .orientation -'111', '110', or '100' crystal
- plane normal aligned with the
- magnetic field, def. is '111'
- .nitrogen -'14N' or '15N', default is '14N'
- sys -Spinach system specification structure
- inter -Spinach interaction specification structure
- Check consistency
