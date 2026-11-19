# test_clexpn_subsystem — native model and retained ownership regression
Keywords: test_clexpn_subsystem cce_compile cce_subsystem units tensors ordering cache configurations meanfield

Registered as kernel/clexpn_subsystem. Tests the native compiler against existing direct algebra and a stock anisotropic Hamiltonian, scalar Hz/rad/s conversion, mixed multiplicities, nontrivial observed order, two observed spins, empty bath clusters, physical/local mappings, isotope/coordinate/group retention, positive nonuniform bath states and common exterior expectations. Cache checks distinguish identical calls from changes to representation, assumptions, model, and bath density. Duplicate IDs, negative populations and stale full-system basis metadata are covered. Returns the standard regression result; no assertion in an existing test is weakened.

## Native interface

Uses native compilation twice for deterministic algebra equality and an independent stock Hamiltonian for physical validation. It does not call an obsolete compiler. All retained physical-order, cache-key, population, and ownership assertions remain.
