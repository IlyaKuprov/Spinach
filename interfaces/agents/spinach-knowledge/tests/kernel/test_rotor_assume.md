# tests/kernel/test_rotor_assume.m

- Signature: `result=test_rotor_assume()`

## Purpose

Explicit assumptions throughout rotor-stack frame construction.

## Physical / mathematical content

An anisotropic heteronuclear pair supplies complex, noncommuting rotor Hamiltonians in both MAS orientation conventions.

## Numerical / algorithmic content

Compares fresh/stale/matching objects, numerical-frame repeatability, inconsistent-request refusal, phase grids, and unchanged ordinary rotor stacks. Covers all already-rotating spin classes under `nmr`, `cavity`, `esr`, `deer`, `deer-zz`, `spin-phonon`, and `qnmr`, including direct `rotframe` rejection and fresh/stale `rotor_stack` requests. Valid nuclear transformations under electron-only rotating sets and quadrupolar nuclear transformations under `qnmr` prevent over-rejection. Complex, noncommuting pseudosecular hyperfine dynamics provide the mixed-frame control. The 164 checks include six exact zero-component controls and 36 direct/fresh/stale numerical-frame refusals for the three carrier-free solid-effect components, for both electron and nuclear targets in Hilbert and Liouville space. The zero-coupling, zero-offset fixture exercises a valid Hermitian component that would otherwise acquire a spurious carrier.

## Syntax

`result=test_rotor_assume()`

## Parameters / inputs

None. The test constructs its own bounded physical fixtures.

## Outputs

`result` is the regression record of checks, messages, and failures; the test runner determines its final status.

## Header notes

The regression is registered in `test_manifest` and is not an optimisation or performance benchmark.
