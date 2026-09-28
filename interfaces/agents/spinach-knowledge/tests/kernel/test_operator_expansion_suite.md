# tests/kernel/test_operator_expansion_suite.m

- Signature: `result=test_operator_expansion_suite()`

## Purpose

Regression test for operator expansion, conversion, allocation, and selected tensor and relaxation helpers. Returns a test result with explanatory messages.

## Checks

- Verifies `hilb2liouv` left and right conversions against Kronecker products and its state-vector conversion against column stacking for a non-diagonal matrix.
- Reconstructs a generic spin-one matrix from `oper2ist` coefficients. Checks that `enlev2ist` selects the bottom spin level or top boson level, and that `ct2ist` and `bos2ist` reconstruct their respective central-transition and creation–annihilation–number operator products.
- Reconstructs a finite oscillator matrix and a population projector from `oper2bm` and `enlev2bm` coefficients. Local reconstruction helpers combine coefficients with the complete `irr_sph_ten` or `boson_mono` basis using zero-based state indices.
- For one spin, checks `unit_oper` and empty `mprealloc` matrix dimensions in Zeeman Hilbert, Zeeman Liouville, and spherical-tensor Liouville formalisms; the last uses the current basis size.
- For two spins in Hilbert space, checks the rank-two, zero-projection `twospinist` tensor against an explicit operator combination. Checks that `lindbladian` gives the requested negative decay rate for a diagonal jump process.