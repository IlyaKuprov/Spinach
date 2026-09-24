# tests/interfaces/test_hfc_isotopes.m

- Source: `tests/interfaces/test_hfc_isotopes.m`
- Signature: `result=test_hfc_isotopes()`
- Total lines: 163

## Purpose

Isotope-resolved hyperfine import from Gaussian and ORCA fixtures.

## Physical / mathematical content

Hyperfine tensors scale with nuclear gyromagnetic ratios, including sign changes and off-diagonal tensor components; same-isotope and NMR imports are controls.

## Numerical / algorithmic content

Exercises source-isotope provenance, conversion before threshold/purge, equality boundaries, partial outputs, and explicit refusal to guess missing provenance.

## Syntax

`result=test_hfc_isotopes()`

## Parameters / inputs

None. The test constructs its own bounded physical fixtures.

## Outputs

`result` is the regression record of checks, messages, and failures; the test runner determines its final status.

## Header notes

The regression is registered in `test_manifest` and is not an optimisation or performance benchmark.
