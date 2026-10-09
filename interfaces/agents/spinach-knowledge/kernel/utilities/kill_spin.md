# kernel/utilities/kill_spin.m

`spin_system=kill_spin(spin_system,hit_list)` removes spins or bosonic modes from the unified particle list. `hit_list` is a vector of valid positive particle indices or a logical mask with one element per particle. The returned system retains the surviving physical model, reindexes particle-dependent interactions and relaxation data, and rebuilds an existing basis rather than attempting to patch compiled descriptors.

## Physical model and chemistry

Particle removal updates isotope identity, labels, multiplicities, magnetogyric ratios, frequencies, coordinates, Zeeman and giant-spin tensors, pair couplings, proximity data, and the isotope hash when present. Per-particle relaxation rates and pair-dependent relaxation parameters are restricted to the survivors; scalar-relaxation source labels are reindexed. Bosonic scalar and pair parameters are restricted consistently, including spin labels nested in coupling and Zeeman modulation derivatives. A system with no remaining bosonic particles loses its mode container; otherwise obsolete mode strengths are cleared.

Chemical substance identities and concentrations are preserved even when a substance loses its last spin. Such a substance becomes a spin-free concentration pool. Reaction matching rows touching a removed spin are deleted, and surviving global spin labels are reindexed. Unmatched retained product spins consequently arrive at identity. Named selector electrons are reindexed, but removing an electron essential to a selector is rejected before mutation. Removing a spin from a substance with user-supplied selector matrices is also rejected: those basis-dependent matrices must first be rebuilt explicitly. `kinetics` recompiles reaction maps after basis reconstruction.

## Basis reconstruction and invalidated assumptions

A compiled basis is rebuilt from retained per-substance input settings. Manual columns and numeric longitudinal/zero-quantum labels follow the surviving particles; an isotope-name filter is removed only when its own substance has no matching spin left. Descriptors, offsets, projectors, and the basis cache hash are regenerated together.

Correlation depths are capped by the remaining local population without increasing already smaller depths. IK-DNP uses surviving electron and nucleus counts; IK-SBS uses mode, mode-plus-spin, and nontrivial-spin counts, matching `basis` validation. If either approximation loses a required particle class, that substance switches to IK-0 at the surviving class depth, at least one, without graph-based connectivity pruning. Other substances keep their own settings. An empty substance uses `none` and retains its unit coordinate.

Connectivity, permutation-symmetry settings, and Hamiltonian assumptions are discarded because they describe the old particle set. Zeeman, coupling, giant-spin, and mode strengths derived from those assumptions are likewise invalidated. Reapply `assume` before constructing a new Hamiltonian; warnings identify cleared assumption information.

Source: [kernel/utilities/kill_spin.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/kill_spin.m). [Wiki](https://spindynamics.org/wiki/index.php?title=kill_spin.m).
