# tests/kernel/test_dynamic_remaining_spectral_suite.m

- Signature: `result=test_dynamic_remaining_spectral_suite()`

## Purpose

Tests remaining spectral, symmetry, and Fokker-Planck utilities. Syntax: result=test_dynamic_remaining_spectral_suite()

## Physical / mathematical content

- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.

## Outputs

- result -regression test result with explanatory messages
- The test checks field-swept eigensystem helpers, rotor-stack assembly,
- permutation-symmetry projectors, and one-dimensional Fokker-Planck
- operators against compact analytical references.

## Implementation structure

- Tests remaining spectral, symmetry, and Fokker-Planck utilities. Syntax:
- result=test_dynamic_remaining_spectral_suite()
- result -regression test result with explanatory messages
- The test checks field-swept eigensystem helpers, rotor-stack assembly,
- permutation-symmetry projectors, and one-dimensional Fokker-Planck
- operators against compact analytical references.
- Announce the test target
- State the utility target of the test
- Check exact diagonalisation path in rspt_eig against eig()
- Check Liouville-space resonance-field extraction for a diagonal pencil
- Check two-root Hilbert-space resonance extraction in a curved level gap
- Check one-dimensional gradient operator construction in Fokker-Planck space
