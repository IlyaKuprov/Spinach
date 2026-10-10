# tests/kernel/test_cwdm_filters.m

Tests the per-substance basis contract using independent three-proton substances and a spin-free pool. Projection, longitudinal, zero-quantum, approximation, manual, and symmetry settings must leave the untouched descriptor identical and retain every unit row. Additional checks cover local manual additions, equivalent isotope/index selectors, IK-DNP depths, invalid input cardinalities, and membership-sensitive hashes.

For a single pair of identical spin-half particles, Hilbert and wavefunction Zeeman symmetry must give triplet/singlet dimensions, while fully symmetric Liouville symmetry has dimension ten. Projector columns are orthonormal to a matrix 1-norm tolerance of 10^-12, and temporary SALC labels leave the compiled Zeeman descriptor cell empty.

A heterogeneous IK-1/none/spin-free fixture checks `space_level={3,[],[]}` against the explicit proximity spelling and verifies that both later descriptors remain unchanged.

A one-spin `create(sys)` call without interaction input must compile one substance at unit concentration and the same basis descriptor as an explicit empty interaction structure.
