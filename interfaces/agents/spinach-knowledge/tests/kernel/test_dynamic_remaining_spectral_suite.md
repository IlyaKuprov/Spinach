# tests/kernel/test_dynamic_remaining_spectral_suite.m

- Signature: `result=test_dynamic_remaining_spectral_suite()`

## Purpose

Tests remaining spectral, symmetry, and Fokker-Planck utilities. Syntax: result=test_dynamic_remaining_spectral_suite()

## Physical / mathematical content

- Tests perturbative energy derivatives and transition properties, resonance-field extraction, one-dimensional gradient and flow operators, rotor phase-grid construction, and S2 permutation-symmetry projector properties.
## Numerical / algorithmic content

- Compares perturbative eigensystem results with exact references; checks one- and two-root resonance fields and Jacobians; compares one-dimensional gradient and flow matrices with references; and tests rotor assembly and symmetry-projector invariants.
## Outputs

- `result` — regression test result with explanatory messages.
## Implementation structure

- Compare `rspt_eig` energies, Hellmann-Feynman derivatives, transition moments, and populations with exact references.
- Check Liouville- and Hilbert-space `eigenfields` resonance extraction, transition moments and widths, population differences, scaled Jacobians, and two-root fields, identities, and Jacobians.
- Compare one-dimensional `g2fplanck` gradient and `v2fplanck` flow operators with matrix references; check `g2fplanck` with empty inactive dimensions.
- Check rank-zero `rotor_stack` assembly and its phase grid.
- Test the S2 fully symmetric projector orthonormality, mixed orbit, and dimension.