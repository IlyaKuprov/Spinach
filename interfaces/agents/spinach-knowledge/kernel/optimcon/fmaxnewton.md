# kernel/optimcon/fmaxnewton.m

- Signature: `[x,data]=fmaxnewton(spin_system,cost_function,guess)`

## Purpose

Finds a local maximum of a function of several variables using Newton and quasi-Newton algorithms. Syntax: [x,data]=fmaxnewton(spin_system,cost_function,guess)

## Physical / mathematical content

- Maximises a user-supplied objective using the method selected in `spin_system.control.method`: LBFGS, regularised BFGS (RBFGS), regularised Newton, or regularised Newton with Goodwin acceleration.
- Search directions use only unfrozen coordinates; after direction construction, the routine applies a line search to select the step.

## Numerical / algorithmic content

- All four methods check only the initial assembled gradient on unfrozen coordinates before constructing a search direction. A norm below `1e-6` is rejected before any Hessian regularisation or solve; the objective value and trajectory shape are not inspected by this guard. This is an optimiser-level check, not a restriction on individual GRAPE contributions.
- Newton and Goodwin request an objective, gradient, and Hessian every iteration. LBFGS and RBFGS request the initial objective and gradient, then reuse line-search gradients. The initial checks use these existing evaluations. With `max_iter=0`, only the objective is evaluated and the initial-guess guard is not applied.
- If the computed direction is not an ascent direction on unfrozen coordinates, the routine falls back to the projected gradient before bracketing and sectioning the line search.

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
