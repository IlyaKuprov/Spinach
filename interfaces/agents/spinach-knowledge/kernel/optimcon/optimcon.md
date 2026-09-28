# kernel/optimcon/optimcon.m

- Signature: `spin_system=optimcon(spin_system,control)`

## Purpose

Validates optimal-control options and updates the spin system object. This call freezes the optimisation problem for ensemble evaluation.

## Parameters / inputs

- `spin_system` — primary Spinach data structure, created by `create.m` and updated by `basis.m` functions.
- `control` — control data structure described in detail in the [online manual](https://spindynamics.org/wiki/index.php?title=optimcon.m). Required fields include `isotopes`, `channels`, `operators`, `rho_init`, `rho_targ`, `pwr_levels`, `pulse_dt`, and `drifts`.

## Outputs

- `spin_system` — updated Spinach data structure.

## Control configuration

- Fidelity measures are `real` (default), `imag`, and `square`. Integrators are `rectangle` (default, piecewise-constant) and `trapezium` (piecewise-linear). Methods are `lbfgs` (default), `rbfgs`, `newton`, and `goodwin`; `newton` and `goodwin` cannot use the trapezium integrator.
- `pulse_dt` specifies positive interval durations. Rectangle integration uses one control value per interval; trapezium integration uses one more value than the number of intervals. Each ensemble drift supplies either one generator or one generator per control value.
- Optional settings include offsets and their operators, power and state ensembles, ensemble correlations and budget, phase cycles, waveform basis or freeze mask, penalties and bounds, distortion functions, Bloch-Siegert corrections, keyholes, trajectory penalties, fidelity timing, plotting, and checkpointing. Unrecognised or mismatched options cause an error.
- `steady` requires `sphten-liouv`; in this mode the supplied `rho_init` is ignored. Hilbert-space and wavefunction controls and drifts must be Hermitian, as must generators used with `goodwin`.
- For `sphten-liouv`, `zeeman-liouv`, and `zeeman-wavef`, nonempty keyhole schedules with `newton` or `goodwin` are not implemented. First-order `lbfgs` and `rbfgs` keyhole methods, empty schedules, and existing Hilbert-space keyhole Hessians remain available. Supplying `control.keyholes` is not supported with trapezium integration; when the field is absent, an empty schedule is created.
- Bloch-Siegert corrections require channel-isotope control operators that are unit quadratures of the canonical `Lx` and `Ly` operators; they are unavailable with `newton`, `goodwin`, or trapezium integration.

## Frozen ensemble problem

`optimcon` builds the ensemble case catalog and records contiguous worker case blocks in `spin_system.control.worker_cases`. With a parallel pool, `SpmdEnabled` must be true. The common frozen problem is published once through `parallel.pool.Constant` in `spin_system.control.invariants`; drift generators are published through `spin_system.control.drift_slices`, using a per-worker `Composite` so each worker receives only the drifts needed for its case block.

Heavy invariants—drift generators, control and offset operators, control commutators, and Bloch-Siegert response operators when present—are removed from the returned control structure; their names are recorded in `spin_system.control.frozen_fields`. The source header documents that other control fields stay live and `ensemble()` re-sends them to workers at every evaluation. When Bloch-Siegert corrections are enabled, channel carrier frequencies are retained in `spin_system.control.carrier_frq`; the header says `bloch_siegert()` uses them for client-side waveform replay. Changes to ensemble composition, operators, generators, channel isotopes, or carrier frequencies require a fresh `optimcon()` call: changing carriers after freezing would make replay disagree with the response operators seen by the optimiser.

## Reference

- [optimcon.m online manual](https://spindynamics.org/wiki/index.php?title=optimcon.m)
