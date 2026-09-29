# kernel/optimcon/grape_liouv.m

- Signature: `[traj_data,fidelity,grad,hess]=grape_liouv(spin_system,drifts,controls,waveform,rho_init,rho_targ,fidelity_type)`
- Source: [kernel/optimcon/grape_liouv.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/grape_liouv.m)

## Purpose

Evaluates the GRAPE objective and its first derivatives for one shaped pulse, propagating an initial state and projecting onto a target state. This is a low-level contribution routine, not the optimiser; its source directs callers to wrappers such as `grape_xy.m` and `grape_phase.m`.

## Inputs and array shapes

- `spin_system` is the Spinach system set up for optimal control. This routine reads settings from `spin_system.control`; it does not establish their defaults.
- `drifts` is a cell array of square drift generators. One matrix represents a time-independent drift; multiple matrices represent time-dependent drifts. `controls` is a cell array of control generators in the selected formalism.
- `waveform` is a real numeric array with one row per control and `spin_system.control.pulse_ntpts` columns; entries are control amplitudes in rad/s. For the rectangle integrator, each column is a piecewise-constant interval. For the trapezium integrator, the columns are pulse nodes and `pulse_dt` has one fewer element than the waveform columns.
- `rho_init` and `rho_targ` are numeric column vectors in the chosen state-vector formalism. The initial state is not used for stroboscopic steady-state optimisation.
- `fidelity_type` selects the real part, imaginary part, or absolute square of the state overlap: `'real'`, `'imag'`, or `'square'`. The supported state-vector formalisms are `sphten-liouv`, `zeeman-liouv`, and `zeeman-wavef`.

The control setting `fid_type` is separate: it must be `'terminal'` or `'average'`. Terminal mode uses the final pulse node; average mode averages the fidelity over nodes 1 through N. The source does not assign a default for either setting.

## Derivatives, masks, and limits

The returned `grad` differentiates the fidelity with respect to the waveform amplitudes. If `spin_system.control.freeze` is empty, the routine substitutes an all-false mask the size of `waveform`; otherwise it indexes the supplied mask by control row and time column. Frozen entries receive zero gradient. When a Hessian is requested, frozen rows and columns are zeroed and their diagonal entries are set to one.

The fourth output `hess` is not available with the piecewise-linear `'trapezium'` integrator, stroboscopic steady states, nonempty keyhole schedules, average-fidelity mode, or trajectory penalties. In the supported state-vector formalisms, keyholes are also rejected for `'newton'` and `'goodwin'` methods; first-order `'lbfgs'` and `'rbfgs'` methods and empty keyhole schedules remain available. These are source guards, not a list of all optimiser methods.

## Additional control transformations

A nonzero `dead_time` pulls the target backward through the last drift. A nonempty `prefix` transforms the initial state using the first drift, and a nonempty `suffix` transforms the target using the last drift. With the Bloch-Siegert option enabled, the per-control response operator contributes an amplitude-squared term to the slice generator, with its corresponding amplitude derivative.

The routine sums the operators in `traj_pen`, averages their expectation values over the same nodes used by the average-fidelity calculation, and subtracts this cost and its gradient. Trajectory penalties are rejected for stroboscopic steady states and whenever `phase_cycle` is nonempty. This function contains no phase-cycle mask or wrapper transformation; those are outside this mapped source. The forward trajectory, when returned, is a stack of state vectors across the pulse nodes. Zero fidelities and derivatives are accepted for auxiliary costates; initial-guess checks belong to the assembled objective rather than this contribution.

[Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=grape_liouv.m)
