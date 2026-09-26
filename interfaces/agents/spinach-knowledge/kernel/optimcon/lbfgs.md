# kernel/optimcon/lbfgs.m

- Signature: `direction=lbfgs(dx_hist,dg_hist,g)`

## Purpose

Calculates an approximation to the Newton-Raphson search direction for maximising a function using past gradients to build a serviceable substitute to a Hessian. The Hes- sian matrix is not formed explicitly. Syntax: direction=lbfgs(dx_hist,dg_hist,g)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.

## Numerical / algorithmic content

## Parameters / inputs

- dx_hist -history of x increments, a stack
- of column vectors, from the latest
- to the earliest
- dg_hist -history of gradient increments,
- a stack of column vectors, from
- the latest to the earliest
- g -current gradient

## Outputs

- direction -LBFGS approximation to the
- maximisation step vector

## Implementation structure

- Calculates an approximation to the Newton-Raphson search
- direction for maximising a function using past gradients
- to build a serviceable substitute to a Hessian. The Hes-
- sian matrix is not formed explicitly. Syntax:
- direction=lbfgs(dx_hist,dg_hist,g)
- dx_hist -history of x increments, a stack
- of column vectors, from the latest
- to the earliest
- dg_hist -history of gradient increments,
- a stack of column vectors, from
- the latest to the earliest
- g -current gradient
