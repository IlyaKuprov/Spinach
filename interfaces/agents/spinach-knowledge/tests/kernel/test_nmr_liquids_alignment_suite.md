# tests/kernel/test_nmr_liquids_alignment_suite.m

- Signature: `result=test_nmr_liquids_alignment_suite()`

## Purpose

Tests compact literature-alignment probes for liquid-state NMR pulse sequences. Syntax: result=test_nmr_liquids_alignment_suite()

## Physical / mathematical content

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Outputs

- result -regression test result with explanatory messages
- The test covers compact gCOSY, HMBC, HSQC, and TOCSY paths that
- were updated during the nmr_liquids literature-alignment pass.

## Implementation structure

- Tests compact literature-alignment probes for liquid-state NMR
- pulse sequences. Syntax:
- result=test_nmr_liquids_alignment_suite()
- result -regression test result with explanatory messages
- The test covers compact gCOSY, HMBC, HSQC, and TOCSY paths that
- were updated during the nmr_liquids literature-alignment pass.
- Announce the test target
- State the test target
- Build a compact homonuclear system for gCOSY
- Set compact gCOSY parameters
- Run the P-type gCOSY pathway
- Run the N-type gCOSY pathway
