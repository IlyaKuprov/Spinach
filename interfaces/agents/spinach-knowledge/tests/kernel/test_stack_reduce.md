# tests/kernel/test_stack_reduce.m

- Signature: `result=test_stack_reduce()`

## Purpose

Regression tests for reduction of horizontal stacks of true state columns, including permutation symmetry, wide sparse inputs, and independently scaled weak columns.

## Physical / mathematical content

A magnetically equivalent proton pair supplies symmetric and antisymmetric sectors in `sphten-liouv`, `zeeman-liouv`, and `zeeman-wavef`. Single, repeated, phase-shifted, mixed-sector, and complex columns are compared with direct matrix exponentiation; Liouville-space generators also include scalar dissipative decay. Projector checks require both state preservation and genuinely smaller state spaces.

## Numerical / algorithmic content

The test compares stacked ZTE with the union of independently screened trajectories and checks zero-generator coordinate masks, explicit `nstates` bounds, and row ranking. Density, small-norm, and explicit-disable shortcuts retain their existing precedence. Disconnected-subspace path tracing must preserve the support of complex columns.

A real process pool is required for the 300-by-40 mixed-scale case, which satisfies the row/column/pool conditions that previously exposed whole-stack propagation scaling. The weak amplitude is `sqrt(zte_tol*eps('double'))`: it is explicitly verified to be above the actual ZTE tolerance and below machine epsilon, rather than assuming epsilon is below the default `1e-24` tolerance. Both sparse and dense stacks must retain the weak column's coupled coordinate. A growing non-unitary mode checks that a nonzero column initially below tolerance is not skipped; a nilpotent generator gives exact dynamical row maxima for `nstates` ranking.

A 4096-by-513 sparse stack, containing 512 phased nonzero columns and one zero column, must reduce to the two reachable rows. This bounded regression exercises streamed propagation without requesting an out-of-memory allocation.

## Syntax

```matlab
result=test_stack_reduce()
```

## Parameters / inputs

None. The test constructs its own systems and requests one process worker through `create`.

## Outputs

`result` is the standard Spinach regression result structure with explanatory messages and a failure list. The registered test identifier is `kernel/stack_reduce`.

## Header notes

Requires the Spinach kernel path and `tests/lib`. The normal test runner supplies the test paths; `run_tests('pattern','kernel/stack_reduce')` selects this regression.
