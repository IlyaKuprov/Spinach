# tests/kernel/test_cwdm_blocks.m

Checks T2 on two independently coupled proton–carbon substances: Hamiltonian,
T1/T2 relaxation, commutation pulse generators, and left/right product actions
have no stored nonzero outside their compiled offsets blocks. Exact and cheap
single-spin detection states have support only in their hosting block. These
are exact sparse-pattern assertions, not tolerance-based norm comparisons.
Identity checks cover individual spins, local products, numeric sums, isotope sums, and all-spin sums, with both `E` and `T0,0` local requests. Exact, cheap, and deprecated chemical identities all carry the hosting substance concentration. These checks do not establish the broader concentration semantics of thermal equilibrium or detection.

Identity operator checks cover left, right, commutator, and anticommutator actions, with individual `E`/`T0,0` requests and local products in CSC and XYZ formats. Isotope, all-spin, and numeric sums check per-spin multiplicity and exact hosting-block support.

Level-projector checks compare both hosting blocks with independently built single-substance `ZL1` states for cheap, exact, and chemical construction, including the identity term and concentration weighting.

Partner-state checks exercise both hosting substances, preserve full global descriptor arrays, compare every returned state with its local product, and reject a genuinely cross-substance partner specification.

Detection support checks use unweighted `coil_state`; identity assertions use concentration-weighted `state` in exact, cheap, and deprecated chem modes.
