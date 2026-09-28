# kernel/optimcon/wrappers/grape_phase.m

- Signature: `[traj_data,fidelity,gradient,hessian]=grape_phase(phi_profile,spin_system)`

Evaluates a GRAPE objective for phase controls with amplitudes fixed by `spin_system.control.amplitudes`. Converts amplitude–phase samples and the freeze mask to Cartesian controls, calls `grape_xy`, then converts requested gradients and Hessians back to phase derivatives. `traj_data` contains trajectory data; `fidelity` is the objective value. Penalties occupy separate slices of fidelity and derivative outputs. The Hessian is reordered for MATLAB vectorisation.

Validation requires real phase and nonnegative real amplitude arrays with equal element counts, half as many rows as the even control count, and columns matching the integrator’s time-step convention.

Source: https://spindynamics.org/wiki/index.php?title=grape_phase.m