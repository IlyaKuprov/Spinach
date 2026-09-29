# kernel/optimcon/grape_hilb.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/grape_hilb.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=grape_hilb.m)

## Purpose and interface

`[traj_data,fidelity,grad,hess]=grape_hilb(spin_system,drifts,controls,waveform,rho_init,rho_targ,fidelity_type)` evaluates a Hilbert-space density-matrix GRAPE objective and requested derivatives. It is a low-level evaluator, not the optimiser; the source directs users to [`grape_xy.m`](wrappers/grape_xy.md) or [`grape_phase.m`](wrappers/grape_phase.md) for the higher-level optimisation interface.

`spin_system` carries settings prepared by `optimcon.m`. `drifts` supplies drift Hamiltonians, cycled across pulse intervals; `controls` is a cell array of Hilbert-space control operators. The real numeric `waveform` has one row per control and `pulse_ntpts` columns, with amplitudes in rad/s. `rho_init` and `rho_targ` are density matrices. The final argument selects `'real'`, `'imag'`, or `'square'` fidelity.

## Propagation and fidelity

Optional dead-time, prefix, and suffix operations are applied around the pulse propagation; per-interval keyholes can also be included. With the `rectangle` integrator, there is one waveform sample and one `pulse_dt` entry per interval. With `trapezium`, the waveform has one more sample than the number of intervals and the generator uses Iserles–Nørsett product quadrature, so an interior sample contributes to both adjacent intervals. Rectangle propagation can include configured Bloch–Siegert response operators, whose first-order contribution is also included in the gradient.

The overlap is evaluated as `hdot(rho_targ,rho_at_node)`, with the target on the conjugated side. The final argument selects its real part, imaginary part, or absolute square. Separately, `spin_system.control.fid_type` selects terminal fidelity or an average over propagated nodes 1 through `nsteps`. A backward costate sweep computes the gradient; for a trapezium waveform, an interior sample contributes through both neighbouring intervals. Zero fidelities and derivatives are valid for auxiliary costates.

If `traj_pen` is configured, the routine sums the penalty operators, evaluates their real expectation over the propagated nodes, and subtracts the mean penalty from fidelity and gradient. Average-node fidelity or trajectory penalties are incompatible with the requested Hessian. Penalties combined with a nonempty `phase_cycle` are rejected. This evaluator does not apply a phase-cycle waveform mask or freeze mask.

## Derivative outputs and constraints

Requesting the third output computes `grad`, with the same shape as `waveform`. The fourth output requests a Hessian, returned as a square matrix of dimension `numel(waveform)`. A Hessian is implemented only for the `rectangle` integrator and requires `spin_system.control.method` to be `'newton'` or `'goodwin'`; the source has no trapezium Hessian path. Hessians with trajectory cost terms are rejected.

The checks require `zeeman-hilb` formalism, square numeric state matrices, dimensionally consistent square drift and control matrices, a real numeric waveform with the required row and time-point counts, valid fidelity selections, and timing-grid lengths consistent with the integrator. The accepted integrators are `'rectangle'` and `'trapezium'`. `traj_data.forward` contains `nsteps+1` forward-state cells only when trajectory return or supported trajectory plotting is requested; otherwise it is empty.
