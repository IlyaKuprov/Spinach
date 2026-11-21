# tests/kernel/test_cwdm_states.m

T6 compares two proton-pair blocks with independent unit-concentration molecules. Exact `coil_state` vectors are identical to embedded single-molecule operator shapes, and `state` multiplies those shapes by the hosting concentration exactly once. The labels include longitudinal and raising operators, identity, and a level projector. The all-spin sum checks separate weighting of both blocks.

`equilibrium` is compared with independently weighted thermal vectors using a 1e-12 absolute vector-norm bound. Zero spin-bearing populations leave coils unchanged, remove their state amplitudes, and leave a populated spin-free substance with only its unit coordinate. The actual stock-layout comparison is external to the registered suite.

Malformed methods, repeated spin indices, mismatched product descriptions, and wrong projection-array lengths must be rejected by the public `state` grumbler.

Single-substance wavefunction construction must return the same unit ket at concentrations zero, 0.3, and two.

The fixed three-argument `coil_state` API is tested: an omitted spin list and a fourth argument are rejected, and an explicit wavefunction call returns the unit ket.

The weighted `state` wrapper rejects the removed shortcut, while `coil_state` rejects a fourth argument, in all four supported formalisms.

Direct exact tensor construction is compared with full operator action on a mixed proton/deuteron/sodium system, in full and single-spin-correlation bases. Products, complex Cartesian operators, higher-rank tensors, central transitions, level projectors, and identities are covered, including descriptors projected out of the truncated basis.
