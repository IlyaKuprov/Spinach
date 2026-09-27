# etc/diamond_defects/diamond_nv0_es.m

- Signature: [sys,inter]=diamond_nv0_es(parameters)

## Purpose

Returns the NV0 excited-state electron–nitrogen spin system for diamond. Magnetic parameters are from Felton et al., *Phys. Rev. B* **77**, 081201 (2008), https://doi.org/10.1103/PhysRevB.77.081201.

## Physical / mathematical content

The effective electron has spin 3/2 (E4), with principal g values 2.0035, 2.0035, and 2.0029 and axial zero-field splitting D = 1685 MHz. The nitrogen hyperfine principal values for 15N are -23.8, -23.8, and -35.7 MHz. For 14N they are scaled by the ratio of nuclear gyromagnetic ratios. No nitrogen quadrupole interaction is included because none is reported for this state.

## Numerical / algorithmic content

The routine constructs the trigonal principal-axis frame, rotates the electron Zeeman, zero-field-splitting, and hyperfine tensors for the requested crystal orientation, and stores them in Spinach's interaction matrices.

## Parameters / inputs

- parameters.orientation: '111', '110', or '100'; the corresponding crystal-plane normal is aligned with the magnetic field.
- parameters.nitrogen: '14N' or '15N'.

## Outputs

- sys: Spinach system specification structure.
- inter: Spinach interaction specification structure.

## Implementation structure

The function validates both input fields, selects the isotope-specific hyperfine tensor, applies the orientation rotation, and returns the E4–nitrogen system with its Zeeman, zero-field-splitting, and hyperfine interactions.
