# kernel/states/coil_state.m

`rho=coil_state(spin_system,states,spins,method)` constructs an unweighted detection operator as a state vector or density matrix. It preserves the state-description API: isotope and numeric spin sums, local product specifications, level projectors, and unweighted per-substance wavefunctions. All four arguments are required: the method is `exact`, and wavefunctions use an empty spin list; `chem` is not an unweighted method and is rejected.

In spherical-tensor Liouville space, exact construction multiplies tensor-descriptor coefficients by the single-spin left-product table entries acting on the unit coordinate, then projects them directly into the hosting basis. It does not allocate full superoperators. Identity terms retain their hosting substance, including in level-projector expansions. Sums contribute once per selected spin; local identity products contribute once. Product specifications crossing substances are rejected by `which_subst`.

Concentrations never enter this constructor: coils for a zero-population substance remain nonzero. Use `state` for concentration-weighted initial states, so coil detection weights a signal once, not twice. Zeeman Liouville construction uses a direct sum of vectorised local identities; Hilbert construction returns the block-local operator.

In wavefunction formalism, this unweighted primitive stacks one product ket per substance, including the scalar ket for a spin-free block. This is storage only: no concentration weighting or incoherent mixing is represented. The concentration-weighted `state` wrapper rejects segmented wavefunction requests.
