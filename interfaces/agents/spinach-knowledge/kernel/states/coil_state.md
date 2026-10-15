# kernel/states/coil_state.m

`rho=coil_state(spin_system,states,spins)` constructs an unweighted detection operator as a state vector or density matrix. It preserves the state-description API: isotope and numeric spin sums, local product specifications, level projectors, and unweighted per-substance wavefunctions. All three arguments are required; wavefunctions use an empty spin list. There is no method selector.

In spherical-tensor Liouville space, exact construction multiplies tensor-descriptor coefficients by the single-spin left-product table entries acting on the unit coordinate, then projects them directly into the hosting basis. It does not allocate full superoperators. Identity terms retain their hosting substance, including in level-projector expansions. Sums contribute once per selected spin; local identity products contribute once. Product specifications crossing substances are rejected by `which_subst`.

Concentrations never enter this constructor: coils for a zero-population substance remain nonzero. Use `state` for concentration-weighted initial states, so coil detection weights a signal once, not twice. Zeeman Liouville construction uses a direct sum of vectorised local identities; Hilbert construction returns the block-local operator.

In wavefunction formalism, this unweighted primitive stacks one product ket per substance, including the scalar ket for a spin-free block. This is storage only: no concentration weighting or incoherent mixing is represented. The concentration-weighted `state` wrapper rejects segmented wavefunction requests.

Direct spherical-tensor projection hashes full query descriptor rows and matches them to cached `bas.row_hash{n}` using `hash_ismember`; reference basis rows are not rehashed.
