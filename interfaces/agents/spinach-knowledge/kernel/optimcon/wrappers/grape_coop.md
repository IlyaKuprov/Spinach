# kernel/optimcon/wrappers/grape_coop.m

- Signature: `[traj_data,fidelity,gradient]=grape_coop(phi_profile,spin_system)`

Builds a cooperative two-pulse phase cycle for a point-to-point transformation. The executable code divides `phi_profile` into row blocks: the first `ncontrols/2` rows feed pulse A and the remaining rows feed pulse B. Each block is evaluated by `grape_phase`; `traj_data` returns the two trajectory outputs as a pair.

For each pulse, the wrapper removes the component of the final state parallel to the single target state. It sums the two residual impurities and subtracts the mean squared norm of that sum from the first fidelity slice; Hilbert-space states use `hdot`, while other formalisms use an explicit orthogonal projector. The reported fidelity averages the two pulse fidelities. The gradient stacks their phase gradients, averages them, and subtracts the impurity-cancellation gradient from the first slice; other penalty slices remain separate.

The guards require exactly one initial state, one nonzero target, an even control count, and reject average-fidelity mode and trajectory penalties. The function is phase-modulated and does not run an optimiser or select a line-search method.

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/wrappers/grape_coop.m) · [Wiki](https://spindynamics.org/wiki/index.php?title=grape_coop.m)
