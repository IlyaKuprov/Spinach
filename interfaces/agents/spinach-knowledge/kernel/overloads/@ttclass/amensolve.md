# kernel/overloads/@ttclass/amensolve.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/amensolve.m`
- Signature: `x=amensolve(A,y,tol,opts,x0)`
- Total lines: 574

## Purpose

Solves the linear system Ax=y using the AMEn iteration. Syntax: x=amensolve(A,y,tol,opts,x0)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`, `revert_interface()`, `reduce_matrix()`, `reduce_vector()`, `local_matvec()`, `local_matrix()`, `local_vector()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- A -ttclass representing a square matrix
- y -ttclass representing the right hand side
- tol -relative approximation and stopping tolerance,
- 1e-6 is a good start
- x0 -ttclass representing the initial guess
- Options (pass empty array for defaults):
- opts.nswp -maximum number of AMEn sweeps
- opts.init_guess_rank -the rank of the initial guess
- opts.enrichment_rank -the rank of the residual and
- enrichment
- opts.resid_damp -local accuracy gap
- opts.rmax -maximum TT rank limit for the solution
- opts.max_full_size -direct vs iterative solver switch-
- over dimension
- opts.local_iters -maximum number of bicgstab iterations
- for local problems
- opts.verb -Verbosity level: silent (0), sweep (1) or full (2)

## Outputs

- x -ttclass representing the solution such that
- |x-X| < tol |X| in Frobenius norm, where X is
- the exact solution.
- Note: A, y and x0 should have ntrains==1. Call shrink()
- on all three if that is not the case.

## Implementation structure

- Solves the linear system Ax=y using the AMEn iteration. Syntax:
- x=amensolve(A,y,tol,opts,x0)
- A -ttclass representing a square matrix
- y -ttclass representing the right hand side
- tol -relative approximation and stopping tolerance,
- 1e-6 is a good start
- x0 -ttclass representing the initial guess
- Options (pass empty array for defaults):
- opts.nswp -maximum number of AMEn sweeps
- opts.init_guess_rank -the rank of the initial guess
- opts.enrichment_rank -the rank of the residual and
- enrichment

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `isfield()`, `exist()`, `ttort()`, `grumble()`, `clearcoeff()`, `local_vector()`, `local_matvec()`, `local_matrix()`, `local_iters()`, `bicgstab()`, `frob_chop()`, `norm_x()`, `reduce_matrix()`, `reduce_vector()`, `norm_yAx()`, `revert_interface()`.
