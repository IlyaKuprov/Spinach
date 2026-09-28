# examples/fundamentals/state_tests/normalization_2.m

- Signature: `normalization_2()`

## Purpose

Compare the squared norms (inner products) of Cartesian spin states across Spinach's Zeeman Hilbert, Zeeman Liouville, and spherical-tensor Liouville formalisms.

## Method

The test uses the same two-spin system as `normalization_1`: `1H` and `235U` at 14.1 T, Zeeman scalars 2.5 and 1.0, and scalar coupling 10. For each formalism it records `trace(A'*A)` for the x, y, and z states of both spins, producing six values per formalism. The Zeeman Hilbert and Zeeman Liouville columns must agree; the spherical-tensor Liouville column is compared after scaling by the product of spin multiplicities. Either comparison fails above `1e-6` (L1 norm). Thus the three columns are not asserted to be identical without that normalization.
