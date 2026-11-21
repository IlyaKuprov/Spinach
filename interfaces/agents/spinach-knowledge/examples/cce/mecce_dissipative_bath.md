# examples/cce/mecce_dissipative_bath.m

Signature: `[coherence,times]=mecce_dissipative_bath()`

Reproduces the five ME-CCE curves in Fig. 1(c) of [Onizhuk et al., Phys. Rev. Lett. 132, 250401 (2024)](https://doi.org/10.1103/PhysRevLett.132.250401). Table S1 in the [APS supplement](https://journals.aps.org/prl/supplemental/10.1103/PhysRevLett.132.250401/supplement.pdf) fixes the probe couplings and the twelve periodic nearest-neighbour XXZ bonds. The bath starts in a Néel state. Equal raising and lowering rates apply to each site: this phenomenological high-temperature relaxation does not derive a microscopic T1. The probe coupling is longitudinal, so Spinach can evolve a projected coherence block rather than the entire probe density. Connected ring clusters retain static fields from the exterior Néel spins.

The returned 201-by-5 array contains normalised coherence magnitudes, with expansion order in columns. The 201-element time row is in the published panel's plotted coordinates: comparing to the published curves fixes the time mapping to the tabulated frequencies used directly, with rate 0.01 for each jump and no division by the realised sample-mean probe coupling. The 101 blue markers are digitised from the published exact calculation and are not a native full-bath simulation; the inset subtracts these samples from each corresponding ME-CCE curve.

The original two-spin, single-raising-jump example setup remains as a test fixture for the dense Liouvillian comparison in `test_cce_examples`.

[MATLAB source](../../../../../examples/cce/mecce_dissipative_bath.m) · [Exact marker fixture](../../../../../examples/cce/mecce_dissipative_bath_exact.csv) · [Independent density regression](../../../../../tests/kernel/test_cce_examples.m)
