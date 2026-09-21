# tests/kernel/test_nmr_liquids_alignment_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_nmr_liquids_alignment_suite.m`
- Signature: `result=test_nmr_liquids_alignment_suite()`
- Total lines: 193

## Purpose

Tests compact literature-alignment probes for liquid-state NMR pulse sequences. Syntax: result=test_nmr_liquids_alignment_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_spin_system()`, `liquid()`, `test_close()`, `all()`, `fid_p()`, `fid_n()`, `isstruct()`, `isfield()`, `fftshift()`, `conj()`, `spec_pn()`, `fid_hmbc()`, `isequal()`, `state()`, `speye()`.
