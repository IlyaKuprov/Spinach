# tests/kernel/test_cce_parallel.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_parallel.m)

## Purpose

Serial, one/two workers, common samples, and resource caps.

## Behaviour

Checks the following named assertions in the source:

- `one-worker raw identity`.
- `one-worker reduction`.
- `two-worker sample identity`.
- `two-worker raw identity`.
- `two-worker reduction`.
- `state scheduling`.
- `analytic sampled signal`.
- `task cap`.
- `memory cap`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_cce_parallel()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status, and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_parallel.m)

Checks resident-memory refusal for 200 observable columns on both gCCE and generalised master-equation paths; the four raw output records alone require 6,020,000 bytes against a 5,000,000-byte cap.

Public optional-field contracts are explicit in reused fixtures: irrelevant method, event, geometry, and legacy numerical options are removed rather than passed silently.

The concurrency fixture closes its own one-worker pool before create requests two workers, and asserts the resulting live worker count. create deliberately reuses an existing pool. ME memory checks use the registered native evaluator with an explicit Liouville representation; the observable-memory test supplies physical bath_jumps.

## Native scheduling and resource validation

The registered regression uses clexpn and its registered experiment evaluator for serial, one-worker and two-worker scheduling. It preserves seeded assignments, fixed-order raw records and reduction, independent analytic Ising phases, concurrent physical-space and Liouville allocation refusals, single-worker admission, multi-observable ledger limits, and retained master-generator storage checks. Projected jump construction is checked through cce_subsystem and cce_generators against the independent decay generator -eye(4)/2. A memory budget that admits one task workspace reduces the worker count and preserves every result; a budget that admits none is refused.

Memory estimates are computed from the actual native model, ensemble, schedule, and parameter structures. Byte counts are schema-dependent bookkeeping, not physical observables; the exact dimensional arithmetic and the accept/refuse decisions remain independently testable.
