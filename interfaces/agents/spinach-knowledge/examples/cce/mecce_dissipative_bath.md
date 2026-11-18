# examples/cce/mecce_dissipative_bath.m

Signature: `[answer,diagnostics]=mecce_dissipative_bath()`

An observed electron and one bath carbon-13 nucleus evolve under the illustrative hyperfine tensor diag([1100 900 2200]) Hz and a nuclear raising-jump dissipator with rate 3000 per second. Standard sys/inter/bas/parameters blocks select a zeeman-liouv retained basis, generalised CCE order 1, transition preparation, density output, and bath_jumps with the physical descriptor 'L+'. clexpn calls free_evolution; the example constructs no state, detector, or jump matrix.

answer is 31-by-4: column-vectorised observed 2-by-2 density matrices, one row per time on 0 to 200 microseconds. The plotted signal is real(answer(:,3)). Uncontrolled generalised Liouville density propagation uses the stock exponential-action step with entry chopping disabled and transports its floating-point estimate. Independent dense Liouvillian/expm propagation and the comparison curve are retained in test_cce_examples. The supplied jump rate is phenomenological, not a microscopically derived T1.

[MATLAB source](../../../../../examples/cce/mecce_dissipative_bath.m) · [Independent reference tests](../../../../../tests/kernel/test_cce_examples.m)
