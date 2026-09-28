# kernel/optimcon/alpha_conds.m

- Signature: `test=alpha_conds(test_type,alpha,fx_0,fx_1,gfx_0,gfx_1,dir,spin_system)`

## Purpose

Applies a selected line-search acceptance test used by the bracketing and sectioning routines in constrained optimisation. Returns true when the selected condition is satisfied.

## Numerical / algorithmic content

- `test_type=0`: monotonic increase, `fx_1 > fx_0`.
- `test_type=1`: Armijo sufficient increase, `fx_1 >= fx_0 + spin_system.control.ls_c1*alpha*(gfx_0'*dir)`.
- `test_type=2`: strong Wolfe curvature, `abs(gfx_1'*dir) <= spin_system.control.ls_c2*abs(gfx_0'*dir)`.
- `test_type=3`: positive directional derivative at the trial point, `gfx_1'*dir > 0`.

## Parameters / inputs

- `test_type` — condition selector (`0`, `1`, `2`, or `3`).
- `alpha` — trial step length; used by test 1.
- `fx_0`, `fx_1` — objective values at the initial and trial points; used by tests 0 and 1.
- `gfx_0` — gradient at the initial point; used by tests 1 and 2.
- `gfx_1` — gradient at the trial point; used by tests 2 and 3.
- `dir` — search direction vector; used by tests 1, 2, and 3.
- `spin_system` — Spinach data structure containing the line-search settings in `control`.

## Outputs

- `test` — logical true if the selected condition is satisfied.

## Implementation structure

The function checks that `test_type` is a real scalar in `0:3`. It validates the inputs used by each test: required objective values and `alpha` are real scalars; required gradients and `dir` are real column vectors with matching dimensions. It then evaluates the selected condition.

[Source documentation](https://spindynamics.org/wiki/index.php?title=alpha_conds.m)