# kernel/optimcon/bfgs.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/bfgs.m`
- Signature: `H=bfgs(dx_hist,dg_hist,g)`
- Total lines: 135

## Purpose

Calculates a BFGS approximation to the Newton-Raphson search direction for maximising a function using past gradients to build a serviceable substitute to a Hessian. Unlike LBFGS, the pseudo-Hessian matrix is formed explicitly. Syntax: H=bfgs(dx_hist,dg_hist,g)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- dx_hist -history of x increments, a stack
- of column vectors, from the latest
- to the earliest
- dg_hist -history of gradient increments,
- a stack of column vectors, from
- the latest to the earliest
- g -current gradient (used for sizing)

## Outputs

- H -BFGS approximation to the Hessian
- matrix corresponding to the *nega-
- tive* Hessian of the objective.
- The corresponding ascent directi-
- on is obtained as: direction=H\g

## Implementation structure

- Calculates a BFGS approximation to the Newton-Raphson search
- direction for maximising a function using past gradients to
- build a serviceable substitute to a Hessian. Unlike LBFGS,
- the pseudo-Hessian matrix is formed explicitly. Syntax:
- H=bfgs(dx_hist,dg_hist,g)
- dx_hist -history of x increments, a stack
- of column vectors, from the latest
- to the earliest
- dg_hist -history of gradient increments,
- a stack of column vectors, from
- the latest to the earliest
- g -current gradient (used for sizing)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `dx_hist()`, `dg_hist()`, `any()`, `g_new()`.
