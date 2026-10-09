# kernel/states/coil_state.m

`rho=coil_state(spin_system,states,spins,method)` constructs an unweighted detection operator as a state vector or density matrix. It preserves the state-description API: isotope and numeric spin sums, local product specifications, level projectors, and single-substance wavefunctions. The optional method is `exact` (default) or `cheap`; `chem` is not an unweighted method and is rejected.

In spherical-tensor Liouville space, an exact state applies the left-product operator to an unweighted unit coordinate; cheap construction locates the requested local tensor descriptor without exact normalisation. Identity terms retain their hosting substance, including in level-projector expansions. Sums contribute once per selected spin; local identity products contribute once. Product specifications crossing substances are rejected by `which_subst`.

Concentrations never enter this constructor: coils for a zero-population substance remain nonzero. Use `state` for concentration-weighted initial states, so coil detection weights a signal once, not twice. Segmented Zeeman construction remains subject to the existing operator/wavefunction rejection boundaries.

Segmented Zeeman-Liouville and wavefunction inputs are rejected before allocation, also when called through the weighted `state` wrapper.
