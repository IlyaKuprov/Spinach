# kernel/optimcon/fmaxnewton.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/fmaxnewton.m`
- Signature: `[x,data]=fmaxnewton(spin_system,cost_function,guess)`
- Total lines: 397

## Purpose

Finds a local maximum of a function of several variables using Newton and quasi-Newton algorithms. Syntax: [x,data]=fmaxnewton(spin_system,cost_function,guess)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.

## Numerical / algorithmic content

- All four methods check only the initial assembled gradient on unfrozen coordinates before constructing a search direction. A norm below `1e-6` is rejected before any Hessian regularisation or solve; the objective value and trajectory shape are not inspected by this guard. This is an optimiser-level check, not a restriction on individual GRAPE contributions.
- Newton and Goodwin request an objective, gradient, and Hessian every iteration. LBFGS and RBFGS request the initial objective and gradient, then reuse line-search gradients. The initial checks use these existing evaluations. With `max_iter=0`, only the objective is evaluated and the initial-guess guard is not applied.
- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `header()`, `footer()`, `itrep()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- spin_system -Spinach data object that has been
- through optimcon.m function
- cost_function -a function handle that takes the input
- the size of guess
- guess -the initial point of the optimisation

## Outputs

- x -the final point of the optimisation
- data.count.iter -iteration counter
- data.count.fx -function evaluation counter
- data.count.gfx -gradient evaluation counter
- data.count.hfx -Hessian evaluation counter
- data.count.rfo -RFO iteration counter
- data.x_shape -output of size(guess)
- data.* -further fields may be set by the
- objective functon
