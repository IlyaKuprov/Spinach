# tests/kernel/test_mode_removal.m

- Signature: `result=test_mode_removal()`

## Purpose

Tests that particle removal preserves retained bosonic interactions and dissipation, resets derived mode assumptions, and permits spin-only rebuilding after the final mode is removed.

## Physical / mathematical content

Hamiltonians and mode dissipators are compared with independently reconstructed retained systems. The fixtures exercise noncommuting mode terms, complex transverse operators, spin-one quadratic terms, pair couplings, first and second spin-modulation derivatives, and mode decay.

## Numerical / algorithmic content

Exercises spin and mode deletion, logical and simultaneous selections, no-op removal, and mode-only retained systems. Compares retained Hamiltonians in exact Zeeman-Hilbert and Zeeman-Liouville formalisms, and finite-temperature dissipators in Liouville space. For C, V, and T modes, checks that removing a spin discards derived mode strengths under `labframe`, `cavity`, and `spin-phonon` assumptions, then compares rebuilt Hamiltonians with independent systems. Final-mode removal from unassumed and previously assumed systems is compared with independently created spin-only systems under `nmr` and `esr` assumptions, using default, `zeeman`, and `couplings` retention in both formalisms.

## Syntax

`result=test_mode_removal()`

## Parameters / inputs

None. The test constructs its own bounded physical fixtures.

## Outputs

`result` is the regression record of checks, messages, and failures; the test runner determines its final status.
