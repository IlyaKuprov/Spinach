# etc/diamond_defects/diamond_r2.m

- Signature: `[sys,inter]=diamond_r2(parameters)`

## Purpose

R2 self-interstitial spin system for diamond. Syntax: [sys,inter]=diamond_r2(parameters) Magnetic parameters from Hunt et al., Phys. Rev. B 61, 3863 (2000), https://doi.org/10.1103/PhysRevB.61.3863

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- parameters is a structure with the following fields:
- .d_sign -sign of D
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- R2 self-interstitial spin system for diamond. Syntax:
- [sys,inter]=diamond_r2(parameters)
- Magnetic parameters from Hunt et al., Phys. Rev. B 61,
- 3863 (2000), https://doi.org/10.1103/PhysRevB.61.3863
- parameters is a structure with the following fields:
- .d_sign -sign of D
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- sys -Spinach system specification structure
- inter -Spinach interaction specification structure
- Check input count
- Check consistency
