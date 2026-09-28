# examples/fundamentals/state_tests/normalization_1.m

- Signature: `normalization_1()`

## Purpose

Check that the generated Cartesian spin states have consistent norms across Spinach's Zeeman Hilbert, Zeeman Liouville, and spherical-tensor Liouville formalisms.

## Method

The test builds a two-spin system (`1H` and `235U`) at 14.1 T, with Zeeman scalars 2.5 and 1.0 and a scalar coupling of 10. In each formalism, it constructs the x, y, and z states for both spins and compares their norms within each set. It fails if any pairwise norm difference exceeds `1e-6`; otherwise it reports that the internal norm consistency test passed.
