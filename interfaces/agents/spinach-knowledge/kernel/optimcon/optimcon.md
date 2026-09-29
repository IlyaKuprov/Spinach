# kernel/optimcon/optimcon.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/optimcon.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=optimcon.m)

## Purpose and call

`spin_system=optimcon(spin_system,control)` validates the optimal-control settings and prepares `spin_system.control` for waveform optimisation. It clears any previous control configuration before processing the new one. It does not evaluate an objective or run an optimisation. The system and basis are produced by `create.m` and `basis.m`; ensemble setup is handled by `ensemble()`, and Bloch–Siegert response operators are built through `bloch_siegert()`.

The routine freezes the ensemble case catalogue and shared invariants for worker evaluation. It places common data and worker-specific drift slices in `parallel.pool.Constant` objects, removes the heavy frozen fields from the live control structure, and records their names in `frozen_fields`. A worker pool, when present, must have `SpmdEnabled`. Control fields not frozen this way remain live for ensemble evaluation. Changes to the ensemble, operators, generators, channel isotopes, or carrier frequencies require a fresh `optimcon` call.

## Required system data and propagation grid

The control structure must provide isotopes, channels, control operators, initial and target states, power levels, pulse timing, and drift generators. The isotope names must be present in the system; channels are positive integer indices into that isotope list, with one channel per control operator. Initial and target states are cell arrays with matching numbers of entries and shapes appropriate to the selected formalism: vectors for Liouville and wavefunction representations, square matrices for `zeeman-hilb`.

`pulse_dt` is a finite, positive row vector. Its sum is `pulse_dur`, and its length is `pulse_nsteps`. The default integrator is `'rectangle'`; it has `pulse_ntpts=pulse_nsteps`. The alternative `'trapezium'` uses `pulse_ntpts=pulse_nsteps+1`. Each ensemble drift entry contains either one generator or `pulse_ntpts` generators, with matrix dimensions matching the control operators. Power levels are required as a finite, positive real row vector.

## Optimisation controls and waveform parameters

`fidelity` selects the overlap projection and defaults to `'real'`; it also accepts `'imag'` or `'square'`. The method defaults to `'lbfgs'`; alternatives are `'rbfgs'`, `'newton'`, and `'goodwin'`. Newton and Goodwin Hessian methods require the rectangle integrator. The iteration defaults are `max_iter=100`, `tol_x=1e-3`, and `tol_g=1e-6`; limited-memory methods default to `n_grads=50` and require at least two history entries. `fmaxnewton` selects bracketing and optional sectioning internally; `optimcon` exposes no line-search selector.

A supplied `basis` must be real numeric, have `pulse_ntpts` columns, and have orthonormal rows within the coded `1e-6` tolerance. The optimiser's guess then contains one coefficient per basis row. A `freeze` mask cannot be combined with a basis; it is converted to logical and checked against the guess by `fmaxnewton`. It may not freeze every waveform coordinate. Without a basis or freeze mask, both fields default to empty.

A supplied `phase_cycle` is a real numeric table, requires an even number of control operators, and has `ncontrols/2+2` columns; its default is `zeros(1,0)`. This is a phase-cycle control table, not the coordinate freeze mask. `optimcon` validates and stores the table but does not transform the waveform. Array-valued bounds are not allowed for phase-modulated optimisation.

## Optional controls and incompatibilities

- `prefix` and `suffix` default to empty and, when supplied, must be function handles. Distortions are function handles in a cell array; absent distortions default to `{@no_dist}`. Distortions are not supported with Newton or Goodwin. `dead_time` defaults to zero.
- `offsets` and `off_ops` must be supplied together as matching cell arrays; otherwise both are empty. Bloch–Siegert corrections default off and are rejected with Newton, Goodwin, or trapezium. When enabled, each control must match the source's unit-quadrature condition for the channel's canonical `Lx`/`Ly` operators, using `1e-6` tests.
- Keyholes default to an empty schedule; a supplied schedule has `pulse_ntpts` entries. Nonempty entries must be function handles, tested for idempotence and linearity at a relative 2-norm tolerance of `1e-6`. The trapezium integrator rejects keyholes. For Liouville and wavefunction formalisms, nonempty keyholes are not implemented with Newton or Goodwin.
- The penalty defaults are `{'SNS'}` with weight `100.0`. Penalties accept `'none'`, `'NS'`, `'SNS'`, `'DNS'`, or `'SNSA'`; weights are non-negative. Bounds default to upper `+1` and lower `-1`. Bounds may be scalar or shaped `[ncontrols,pulse_ntpts]`; array bounds are rejected with `'SNSA'` penalties or phase-modulated optimisation.
- `fid_type` defaults to `'terminal'`; `'average'` is incompatible with Newton/Goodwin and stroboscopic steady states. Trajectory penalties are incompatible with Newton/Goodwin, stroboscopic steady states, and phase cycling. `steady` defaults to false and is available only in `sphten-liouv` formalism. `fidelity` and `fid_type` are distinct: the former chooses the overlap projection, the latter terminal or average node selection.

`ens_corrs` accepts `rho_ens`, `rho_drift`, or `power_drift` and defaults to empty. `rho_ens` cannot be combined with another entry; `rho_drift` requires the number of initial states to equal `ndrifts`, and `power_drift` requires the number of power levels to equal `ndrifts`. `budget` defaults to `Inf`; a finite budget above one must be integral.

`traj_pen` defaults to empty; if supplied, it is a nonempty cell array of state-sized numeric column vectors for Liouville representations or square numeric matrices matching the state dimension for Hilbert/wavefunction representations. The latter are checked for Hermiticity. Plot settings and `distplot` default to empty; `distplot` entries are function handles, and `traj_opts` accepts only `average`. `video_file` is checked when plotting is enabled, and `checkpoint` must be a character string. If `amplitudes` is supplied, it is stored here and validated in `grape_phase`. `parameters` is accepted without parsing; unrecognised remaining control fields cause an error.

The higher-level wrappers [`grape_xy.m`](wrappers/grape_xy.md) and [`grape_phase.m`](wrappers/grape_phase.md) provide pulse parameterisations around the low-level evaluator [`grape_hilb.m`](grape_hilb.md). `optimcon` validates their settings; it does not itself perform waveform optimisation.
