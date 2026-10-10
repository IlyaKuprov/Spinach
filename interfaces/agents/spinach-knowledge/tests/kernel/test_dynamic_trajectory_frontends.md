# tests/kernel/test_dynamic_trajectory_frontends.m

**Source**: [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_trajectory_frontends.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_trajectory_frontends.m)

## Purpose

Regression test for the dynamic trajectory-analysis front-end kernels. It exercises the plotting branches of `trajan()` and the scoring branches of `trajsimil()` on a compact two-spin spherical-tensor trajectory, checking that both helpers expose deterministic branch outputs.

## What the suite checks

- On a two-spin trajectory, `trajan` must plot separate one- and two-spin correlation orders; when an explicit three-point time axis is given, the plotted x-data must match it. Its coherence-order branch must produce five orders from −2 through +2, its per-spin population branch two traces, and its Zeeman level-population branch four traces.
- For identical trajectories, `trajsimil` in RSP mode must return each time point's complex scalar product. RDN, sign-grouped RDN and broad-state-grouped RDN must all return unit similarity at each point. These checks compare returned scores rather than asserting a particular appearance for the plots.
- Plotting tests use invisible figures and restore the previous figure-visibility default; the suite reports failures through the regression result rather than claiming a measured physical trajectory.

## Inputs and outputs

**Syntax**:

```matlab
result = test_dynamic_trajectory_frontends()
```

**Outputs**:

- `result` — regression test result with explanatory messages, accumulated from `test_true` and `test_close` assertions across the `trajan()` and `trajsimil()` branches.

The function takes no inputs.

## References

- `trajan` — trajectory plotting front end exercised in correlation-order, coherence-order, total-each-spin, local-each-spin, and level-population modes.
- `trajsimil` — trajectory similarity scoring front end exercised in RSP, RDN, SG-RDN, and BSG-RDN modes.
- `new_test_result`, `test_true`, `test_close` — regression test harness helpers.
- `test_spin_system`, `unit_state`, `state` — spin system and state construction helpers used to build the test trajectory.
