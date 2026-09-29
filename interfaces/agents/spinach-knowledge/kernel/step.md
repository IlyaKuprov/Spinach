# kernel/step.m

- Signature: `rho=step(spin_system,L,rho,time_step)`

## Purpose

Applies the action of an exponential for one time step without generally forming the full exponential. `time_step` is used directly in products with `L`; this source does not state an independent physical unit for it.

## Generator and state forms

A single `L` gives the one-point, piecewise-constant rule. A two-matrix cell `{L_left,L_right}` uses the piecewise-linear rule; a three-matrix cell `{L_left,L_mid,L_right}` uses the piecewise-quadratic rule. For a manually assembled generator, the source documents `L=H+1i*R+1i*K`, where `H`, `R`, and `K` are the Hamiltonian commutation, relaxation, and kinetics superoperators. A three-element cell whose first entry is a function handle is passed to [`iserstep.m`](pulses/iserstep.md) for state-dependent evolution; the source header identifies `L{2}` as the current time and `L{3}` as the `iserstep` method.

The vector-action form is used for `sphten-liouv`, `zeeman-liouv`, and `zeeman-wavef`; other formalisms use left and right action on a density matrix. For small matrices (below `spin_system.tols.small_matrix`), the code forms `expm(-1i*L*time_step)` and applies it directly: a vector receives left action, while a density matrix receives `P*rho*P'`. Larger cases are subdivided using `ceil(cheap_norm(L)*abs(time_step)/2)` steps and processed by the Taylor-times-vector or commutator-series routine. Numeric states or cells of states are accepted; stack processing may use `parfor`.

## Reporting and side effects

The function normally returns the propagated state without a numeric progress summary. It calls `report` with an all-zero-state warning and returns early for a zero state; when more than 100 substeps are required it reports the substep count and suggests `evolution()`. More than `1e4` substeps raises an error. The Taylor-times-vector routine also calls MATLAB `warning` if its Taylor index exceeds 32. These messages carry counts but no physical units.

When GPU support is enabled in `spin_system.sys.enable`, the function converts numeric generators and states to GPU arrays as needed. This direct implementation does not write files or assign into `spin_system`; it may use GPU arrays and a parallel pool for computation. The state-dependent branch delegates to `iserstep.m`.

## Links

- Related solvers: [`iserstep.m`](pulses/iserstep.md) and [`evolution.m`](evolution.md).
- Source: [`kernel/step.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/step.m).
- Wiki: [`step.m`](https://spindynamics.org/wiki/index.php?title=step.m).
