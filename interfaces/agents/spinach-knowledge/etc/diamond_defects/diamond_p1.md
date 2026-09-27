# etc/diamond_defects/diamond_p1.m

- Signature: [sys,inter]=diamond_p1(parameters)

## Purpose

Returns the P1 centre spin system in diamond. Magnetic parameters are taken from Nir-Arad et al., *Phys. Chem. Chem. Phys.* **26**, 27633 (2024), https://doi.org/10.1039/d4cp03055a, and Smith et al., *Phys. Rev.* **115**, 1546 (1959), https://doi.org/10.1103/PhysRev.115.1546.

## Physical / mathematical content

The electron is coupled to the nitrogen nucleus. For 14N, the principal hyperfine values are 81.3, 81.3, and 114.0 MHz, with a nitrogen quadrupole tensor having D = -3.97 MHz. For 15N, the hyperfine values are -114.0, -114.0, and -159.9 MHz, and no quadrupole tensor is included. The electron g values are 2.00220, 2.00220, and 2.00218.

## Numerical / algorithmic content

The tensors are defined in the trigonal frame and rotated to the chosen crystal orientation. When omitted, orientation and isotope default to '111' and '14N'.

## Parameters / inputs

- parameters.orientation: '111', '110', or '100'; default '111'.
- parameters.nitrogen: '14N' or '15N'; default '14N'.

## Outputs

- sys: Spinach system specification structure.
- inter: Spinach interaction specification structure.

## Implementation structure

The routine validates the structure and applies the isotope-specific electron–nitrogen coupling and orientation rotation before constructing the Spinach interaction matrices.
