# tests/kernel/test_dynamic_remaining_spectral_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_remaining_spectral_suite.m`
- Signature: `result=test_dynamic_remaining_spectral_suite()`
- Total lines: 226

## Purpose

Tests remaining spectral, symmetry, and Fokker-Planck utilities. Syntax: result=test_dynamic_remaining_spectral_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- The file also defines local helper function(s): `local_created_system()`, `local_minimal_system()`, `local_tiny_rank_one()`, `gcp()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `local_minimal_system()`, `rspt_eig()`, `V_ref()`, `test_close()`, `local_tiny_rank_one()`, `orientation()`, `eigenfields()`, `spin()`, `test_true()`, `isequal()`, `local_created_system()`, `g2fplanck()`, `hamiltonian()`, `assume()`, `spdiags()`.
