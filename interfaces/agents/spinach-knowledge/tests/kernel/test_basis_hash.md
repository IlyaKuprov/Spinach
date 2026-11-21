# tests/kernel/test_basis_hash.m

`test_basis_hash()` checks stored row digests against direct MD5 hashes of full local descriptors, dense-query membership, duplicate and missing matches, empty arrays, and spin-free zero-column rows. A mixed spin-half/spin-one multi-substance system supplies full and IK-0-truncated bases. Operators on the truncated basis are compared with full-basis projections for left, right, commutator, and anticommutator actions, including product operators; exact detection states and a truncated reaction-row map are checked separately. Rebuilding as a Zeeman Hilbert context and converting with `sim2liouv` verifies that obsolete spherical-tensor row hashes do not survive.

Run with `run_test('kernel/basis_hash')`.
