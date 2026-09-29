# kernel/optimcon/wrappers/grape_phase.m

- Signature: `[traj_data,fidelity,gradient,hessian]=grape_phase(phi_profile,spin_system)`

Evaluates phase-only GRAPE derivatives with amplitudes fixed by `spin_system.control.amplitudes`. Each phase row pairs with one amplitude row; polar samples are converted to Cartesian x/y controls before calling `grape_xy`. The freeze mask is expanded to both Cartesian channels for each phase coordinate. This is a coordinate wrapper, not an optimiser or line-search selector.

The amplitude and phase arrays must be real numeric arrays with equal element counts and `ncontrols/2` rows; amplitudes must be nonnegative. Under the `rectangle` integrator, each has `pulse_nsteps` columns; under `trapezium`, each has `pulse_nsteps+1`. The control amplitude profile must exist, and `ncontrols` must be even.

Two outputs request trajectory and fidelity only; three also returns the phase gradient; four request gradient and Hessian. The gradient is converted from Cartesian derivatives with `cartesian2polar` and has one third-dimension slice for the objective and each penalty. The Hessian has shape `numel(phi_profile)` by `numel(phi_profile)` by (number of penalties+1): Cartesian channel blocks are reordered, contracted into phase-pair blocks, and the same-track diagonal includes the second-derivative correction from the polar map. Final ordering follows MATLAB vectorisation. The source rejects unsupported output counts and unknown time-propagation algorithms.

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/wrappers/grape_phase.m) · [Wiki](https://spindynamics.org/wiki/index.php?title=grape_phase.m)
