# kernel/optimcon/wrappers/grape_coop.m

- Signature: `[traj_data,fidelity,gradient]=grape_coop(phi_profile,spin_system)`

Evaluates two phase-modulated pulses whose orthogonal target-state impurities are intended to cancel. `phi_profile` contains their phase profiles stacked by control channel. The wrapper runs `grape_phase` for each pulse, returns both trajectories, and averages their fidelities and gradients; it subtracts the mean squared norm of the summed impurities and its gradient from the fidelity slice.

Requires one initial state, one nonzero target, and an even number of controls. Average fidelity type and trajectory penalties are unsupported.

[Source](https://spindynamics.org/wiki/index.php?title=grape_coop.m)