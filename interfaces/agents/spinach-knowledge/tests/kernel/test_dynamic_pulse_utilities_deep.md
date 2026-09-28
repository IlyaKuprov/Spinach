# tests/kernel/test_dynamic_pulse_utilities_deep.m

- Signature: `result=test_dynamic_pulse_utilities_deep()`

Regression checks for pulse utilities, including one-proton Liouville-space cases. The suite compares `grad_pulse` and `grad_sandw` with small-matrix exponential references; checks heterodyne output and RLC responses, including `pwl_tsc`/`pwl` equivalence; inspects point count, shape length, and terminator in a temporary Bruker file; compares finite-RF R-sequence propagators and index maps with explicit references; and checks waveform-basis dimensions and orthonormality and pulse shapes against formulas. It returns a result with explanatory messages.
