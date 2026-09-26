# examples/nmr_solids/case_studies/mathies_carbonate/cp_mas_powder_mhc_fplanck_exchange.m

- Signature: `cp_mas_powder_mhc_fplanck_exchange()`

## Purpose

Cross-polarisation contact curve under magic angle spinning in the presence of chemical exchange for H1, H4 and C19 in the unit cell of monohydrocalcite. Further details in: Calculation time: hours, much faster on a GPU.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Cross-polarisation contact curve under magic angle spinning
- in the presence of chemical exchange for H1, H4 and C19 in
- the unit cell of monohydrocalcite. Further details in:
- Calculation time: hours, much faster on a GPU.
- 400 MHz NMR
- Read CASTEP file
- Drop O and Ca atoms
- Two chemical endpoints: H1, H4, and C19,
- with H1 and H4 under chemical exchange
- Convert shielding tensors into shift using the
- parametrisation of Huang et al. ACIE 2021
- Cartesian coordinates
