# tests/interfaces/test_orca_parser.m

- Signature: `result=test_orca_parser()`

## Purpose

Tests the ORCA log parser using the logs bundled with the examples.

## Physical / mathematical content

For the methyl-radical vacuum DFT example, checks the parsed `g`-tensor and hyperfine tensors against the values ORCA prints in the same log.

## Numerical / algorithmic content

Checks ORCA version detection and maps ORCA's zero-based and incompletely printed nucleus numbering to the corresponding atoms in the coordinate table.

## Parameters / inputs

The test reads ORCA example logs bundled with Spinach.

## Outputs

`result` contains the regression test result and explanatory messages.
