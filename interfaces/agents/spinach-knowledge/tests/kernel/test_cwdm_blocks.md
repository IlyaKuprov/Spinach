# tests/kernel/test_cwdm_blocks.m

Checks T2 on two independently coupled proton–carbon substances: Hamiltonian,
T1/T2 relaxation, commutation pulse generators, and left/right product actions
have no stored nonzero outside their compiled offsets blocks. Exact and cheap
single-spin detection states have support only in their hosting block. These
are exact sparse-pattern assertions, not tolerance-based norm comparisons.
Concentration weighting is outside this test's scope.
