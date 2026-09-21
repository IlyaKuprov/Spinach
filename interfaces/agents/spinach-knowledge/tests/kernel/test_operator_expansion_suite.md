# tests/kernel/test_operator_expansion_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_operator_expansion_suite.m`
- Signature: `result=test_operator_expansion_suite()`
- Total lines: 156

## Purpose

Tests operator expansion and conversion helpers. Syntax: result=test_operator_expansion_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file also defines local helper function(s): `ist_reconstruct()`, `bm_reconstruct()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Outputs

- result -regression test result with explanatory messages
- The test checks that irreducible spherical tensor and bosonic monomial
- expansion helpers reconstruct explicit matrices, and that operator-sized
- allocation helpers return the correct formalism dimensions.

## Implementation structure

- Tests operator expansion and conversion helpers. Syntax:
- result=test_operator_expansion_suite()
- result -regression test result with explanatory messages
- The test checks that irreducible spherical tensor and bosonic monomial
- expansion helpers reconstruct explicit matrices, and that operator-sized
- allocation helpers return the correct formalism dimensions.
- Announce the test target
- State the expansion target of the test
- Check Hilbert-to-Liouville vectorisation identities on a non-diagonal matrix
- Check IST expansion of a generic spin-one matrix
- Check spin and boson energy-level counting conventions in IST expansions
- Check central-transition and boson-product IST expansion wrappers

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `speye()`, `test_close()`, `hilb2liouv()`, `transpose()`, `oper2ist()`, `ist_reconstruct()`, `enlev2ist()`, `ct2ist()`, `centrans()`, `centran()`, `weyl()`, `bos2ist()`, `oper2bm()`, `bm_reconstruct()`, `enlev2bm()`.
