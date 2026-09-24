# kernel/optimcon/optimcon.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/optimcon.m`
- Signature: `spin_system=optimcon(spin_system,control)`
- Total lines: 1619

## Purpose

Validates optimal control options and updates the spin system object. Syntax: spin_system=optimcon(spin_system,control)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.

## Numerical / algorithmic content

- Nonempty keyhole schedules with `newton` or `goodwin` are explicitly not implemented in `sphten-liouv`, `zeeman-liouv`, or `zeeman-wavef`. First-order `lbfgs`/`rbfgs` keyhole methods, empty schedules, and existing Hilbert-space keyhole Hessians remain available; no algorithm is substituted. The same method restriction is enforced at direct `grape_liouv` entry regardless of the requested output count.

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `check_hermiticity()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- spin_system -primary Spinach data structure,
- created by create.m and updated
- by basis.m functions
- control -control data structure described
- in detail in the online manual

## Outputs

- spin_system -updated Spinach data structure
- Note: this function freezes the optimisation problem. The ensemble
- case catalog is built here, its cases are assigned to the
- parallel pool workers in contiguous blocks recorded in
- spin_system.control.worker_cases, and the frozen problem is
- published to the workers exactly once: the common part as a
- parallel.pool.Constant in spin_system.control.invariants, the
- drift generators as a parallel.pool.Constant built from a per-
- worker Composite in spin_system.control.drift_slices, so that
- each worker receives the drifts of its own case block and no-
- thing else; the pool must therefore have SpmdEnabled set to
- true. Heavy invariants -the drift generators, the control
- operators, the offset operators, the control commutators,
- and the Bloch-Siegert response operators -are then removed
- from the returned structure, and their names are recorded
- in spin_system.control.frozen_fields. All other control
- fields stay live: ensemble() re-sends them to the workers
- at every evaluation, and they may be overwritten between
- optimisations. Among them are the Bloch-Siegert channel
- carrier frequencies, kept in spin_system.control.carrier_-
- frq, from which bloch_siegert() rebuilds the response ope-
- rators when a waveform is replayed on the client. Changes
- to the ensemble composition, the operators, the generators,
- the channel isotopes, or the carrier frequencies require a
- fresh optimcon() call: the frozen response operators are
- built from the carriers seen here, and editing them after-
- wards would replay physics that the optimiser never saw.

## Implementation structure

- Validates optimal control options and updates the spin system
- object. Syntax:
- spin_system=optimcon(spin_system,control)
- spin_system -primary Spinach data structure,
- created by create.m and updated
- by basis.m functions
- control -control data structure described
- in detail in the online manual
- spin_system -updated Spinach data structure
- Note: this function freezes the optimisation problem. The ensemble
- case catalog is built here, its cases are assigned to the
- parallel pool workers in contiguous blocks recorded in

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isfield()`, `rmfield()`, `report()`, `banner()`, `ischar()`, `pad()`, `rho()`, `ismember()`, `strcmp()`, `iscell()`, `all()`, `cellfun()`, `isvector()`, `any()`, `check_hermiticity()`.
