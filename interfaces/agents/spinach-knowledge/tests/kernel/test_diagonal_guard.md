# tests/kernel/test_diagonal_guard.m

- Source: `tests/kernel/test_diagonal_guard.m`
- Signature: `result=test_diagonal_guard()`
- Total lines: 144

## Purpose

Supported formalism boundaries for diagonal relaxation retention.

## Physical / mathematical content

Trace preservation, identity stationarity, Hermiticity, and longitudinal/transverse decay rates are checked for spin-half and spin-one baths. A six-proton subsystem from the shipped `molecule_b.xml` additionally checks the imported pure-damping generator in both full Liouville bases, including the Lorentzian relation `R2=pi*FWHM` with FWHM in hertz.

## Numerical / algorithmic content

Requires explicit rejection of Zeeman diagonal retention while retaining spherical diagonal, Zeeman full, and omitted-retention controls. GISSMO full retention must preserve the formerly supported spherical diagonal-retention result exactly.

## Syntax

`result=test_diagonal_guard()`

## Parameters / inputs

None. The test constructs bounded bath fixtures and imports subsystem 2 from `examples/nmr_metabol/molecule_b.xml`.

## Outputs

`result` is the regression record of checks, messages, and failures; the test runner determines its final status.

## Header notes

The regression is registered in `test_manifest` and is not an optimisation or performance benchmark.
