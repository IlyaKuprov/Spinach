# tests/interfaces/test_orca_parser.m

- Signature: `result=test_orca_parser()`

## Purpose

Tests the ORCA log parser on the logs bundled with the examples. Syntax: result=test_orca_parser()

## Physical / mathematical content

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test checks version detection, the mapping of ORCA's zero-based
- and incompletely printed nucleus numbering onto the atoms of the
- coordinate table, and the numerical values of the g-tensor and of
- the hyperfine tensors against the quantities that ORCA itself
- prints separately in the same logs.

## Implementation structure

- Tests the ORCA log parser on the logs bundled with the examples. Syntax:
- result=test_orca_parser()
- result -regression test result with explanatory messages
- The test checks version detection, the mapping of ORCA's zero-based
- and incompletely printed nucleus numbering onto the atoms of the
- coordinate table, and the numerical values of the g-tensor and of
- the hyperfine tensors against the quantities that ORCA itself
- prints separately in the same logs.
- Announce the test target
- State the physical target of the test
- Locate the example logs bundled with Spinach
- Methyl radical, a vacuum DFT calculation with a g-tensor and hyperfines
