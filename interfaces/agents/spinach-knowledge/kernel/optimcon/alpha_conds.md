# kernel/optimcon/alpha_conds.m

- Signature: `test=alpha_conds(test_type,alpha,fx_0,fx_1,gfx_0,gfx_1,dir,spin_system)`

## Purpose

Applies one of the line search acceptance tests used by the brac- keting and sectioning routines in constrained optimisation and returns true when the chosen condition is satisfied. Syntax: test=alpha_conds(test_type,alpha,fx_0,fx_1,... gfx_0,gfx_1,dir,spin_system)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.

## Numerical / algorithmic content

## Parameters / inputs

- test_type -condition selector:
- 0 for monotonic increase test
- 1 for Armijo sufficient increase test
- 2 for strong Wolfe curvature test
- 3 for ascent direction test
- alpha -trial step length
- fx_0 -objective value at the initial point
- fx_1 -objective value at the trial point
- gfx_0 -gradient at the initial point
- gfx_1 -gradient at the trial point
- dir -search direction vector
- spin_system -Spinach data structure with line
- search settings in control

## Outputs

- test -logical true if the selected
- condition is satisfied

## Implementation structure

- Applies one of the line search acceptance tests used by the brac-
- keting and sectioning routines in constrained optimisation and
- returns true when the chosen condition is satisfied. Syntax:
- test=alpha_conds(test_type,alpha,fx_0,fx_1,...
- gfx_0,gfx_1,dir,spin_system)
- test_type -condition selector:
- 0 for monotonic increase test
- 1 for Armijo sufficient increase test
- 2 for strong Wolfe curvature test
- 3 for ascent direction test
- alpha -trial step length
- fx_0 -objective value at the initial point
