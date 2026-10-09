# tests/kernel/test_cwdm_selectors.m

Tests Haberkorn and Jones–Hore selector loss against independent Hilbert-space projectors on a generic positive radical-pair density matrix. It checks nuclear product arrival after electronic partial tracing, channel-resolved concentration balance, and equivalence of named selectors to user-supplied local left/right product matrices.

The integrated CIDNP transport test adds a noncommuting electron–nuclear Hamiltonian. An independent augmented Hilbert generator integrates the singlet-projected nuclear source and the triplet sink; the direct-sum trajectory and nuclear product observable must agree to `1e-10` in the test norm. Converter blocks are divided by their Hilbert dimensions explicitly, so physical traces equal concentrations.

This is a small analytic T14/T15 acceptance case, not the seven-example WP0 radical-pair regression.

The same selective product-arrival trajectory is checked through `evolution` with default reductions and with ZTE enabled. Initially empty product blocks must remain reachable; the comparison detects inappropriate factorisation into independent substance irreps.

Detection and reference operator vectors explicitly use the `exact` method of the four-argument `coil_state` primitive.

The RYDMR boundary checks reject tracked products by identifier and retain the analytic unit singlet yield for untracked recombination without spin mixing.
