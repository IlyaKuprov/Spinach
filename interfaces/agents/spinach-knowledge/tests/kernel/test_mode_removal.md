# tests/kernel/test_mode_removal.m

- Source: `tests/kernel/test_mode_removal.m`
- Signature: `result=test_mode_removal()`
- Total lines: 229

## Purpose

Retained bosonic interactions, dissipation, and assumption lifecycle after particle removal.

## Physical / mathematical content

Hamiltonians and mode dissipators must match independently reconstructed retained systems, including noncommuting and complex couplings.

## Numerical / algorithmic content

Exercises spin/mode deletion, logical and simultaneous selections, no-op removal, mode-only retained systems, pair channels, and nested first/second spin-modulation derivatives. Checks that labframe, cavity, and spin-phonon strengths are discarded and rebuilt for C, V, and T particles. Final-mode removal from unassumed and previously assumed systems is compared with independently created spin-only systems under nmr/esr assumptions and both retention options in exact Hilbert and Zeeman-Liouville formalisms.

## Syntax

`result=test_mode_removal()`

## Parameters / inputs

None. The test constructs its own bounded physical fixtures.

## Outputs

`result` is the regression record of checks, messages, and failures; the test runner determines its final status.

## Header notes

The regression is registered in `test_manifest` and is not an optimisation or performance benchmark.
