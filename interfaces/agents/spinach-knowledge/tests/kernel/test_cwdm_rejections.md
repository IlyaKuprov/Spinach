# tests/kernel/test_cwdm_rejections.m

Checks T3: a scalar coupling between substances raises
`Spinach:create:crossSubstanceCoupling`; product operators and states spanning
substances raise `Spinach:which_subst:crossSubstance`, including an explicit
identity factor. Both identifiers and messages are asserted. An out-of-range
local symmetry index is rejected with the local spin-label error message.

Both Hilbert and Liouville segmented Zeeman systems must raise the named `segmentedZeeman` error from `correlation`, `decouple`, and `homospoil`, as well as from `basis` when permutation symmetry is requested.

Single-substance units are compared exactly with their stock normalisations; single-substance equilibrium is compared with an explicit trace-normalised Boltzmann exponential at absolute and relative whole-array tolerances of `1e-10`.

A two-substance wavefunction fixture with local dimensions two and four requires `Spinach:state:segmentedZeeman`; the corresponding single-substance fully polarised three-spin wavefunction is compared exactly with its eight-component reference.

A disjoint but incomplete chemical partition must raise `Spinach:basis:incompletePartition` rather than silently omit the unassigned spin.

Complete compiled column-vector partitions must yield the same descriptors as row-vector partitions.

Coherent-state tests use local Hilbert dimensions two and three in both Zeeman formalisms, asserting the named segmented rejection and comparing the corresponding single-substance state with the explicit normalised truncated coherent product.

Mixed left-product/commutator equilibrium inputs are rejected in both block orderings with `Spinach:equilibrium:notLeftProduct` and the invalid substance number. Valid left-product blocks are compared against independently normalised Boltzmann exponentials at absolute and relative whole-vector tolerances of `1e-14`.

Exactly zero spinful Hamiltonian blocks in both positions retain the local unit state while the other block retains its independently normalised Boltzmann state.

Cross-substance Hamiltonian entries in either direction, supplied through `I` or the oriented `Q`, must raise `Spinach:equilibrium:crossSubstanceHamiltonian`.

Cross-substance generator entries in either direction must raise `Spinach:reduce:crossSubstanceGenerator` through both `reduce` and `evolution`. With numerical clean-up disabled, supported block-diagonal evolution is compared with the full matrix exponential at absolute and relative whole-vector tolerances of `1e-12`.
