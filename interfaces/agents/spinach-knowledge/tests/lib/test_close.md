# tests/lib/test_close.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_close.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_close.m)

## Purpose

`test_close` adds a numerical regression check with tolerances and an explanation to a test result structure. It compares a value produced by Spinach against a reference value supplied by the test, using an absolute and a relative tolerance, and records the outcome of the check in the result structure.

## Behaviour

- Syntax: `result=test_close(result,label,observed,reference,abs_tol,rel_tol,why)`.
- A consistency-enforcement subfunction (`grumble`) validates all inputs and calls `error` on invalid arguments: `result` must be a scalar structure with `messages` and `failures` fields; `label` and `why` must be non-empty character row strings; `observed` and `reference` must be logical or numeric; `abs_tol` and `rel_tol` must be non-negative real scalars.
- Sparse inputs are converted to full with `full`, and both `observed` and `reference` are converted to double column vectors via `double(...(:))` for the numerical comparison.
- Conditions that prevent the comparison from being made at all are recorded as failures in the same way as a tolerance violation: a size mismatch between `observed` and `reference`, NaN or Inf in the observed vector, or NaN or Inf in the reference vector.
- Where the comparison is possible, the error norm is `norm(observed_vec-reference_vec,2)` and the reference norm is `max([1 norm(reference_vec,2)])`; the tolerance limit is `abs_tol+rel_tol*ref_norm`.
- If the error norm or the limit is non-finite, the check is recorded as a failure with the detail `comparison produced a non-finite scalar`.
- If `error_norm>limit`, the check fails and the failure detail includes the label, the error, the limit, and the `why` explanation.
- A failed check is recorded in the `messages` and `failures` fields of the result structure rather than thrown, so that the checks that follow it in the calling test are still evaluated and the earlier passes are not lost; `run_tests` inspects the `failures` field to decide the status.
- On success, a `PASS` message is appended to `result.messages` containing the label, the error norm, the tolerance limit, and the `why` explanation; on failure, a `FAIL` message is appended to `result.messages` and the detail is appended to `result.failures`.
- The function calls `test_record(result)` to retain the record for the `run_tests` catch path.

## Inputs and outputs

**Inputs**

- `result` — test result structure.
- `label` — check label.
- `observed` — value produced by Spinach.
- `reference` — comparison value supplied by the calling test.
- `abs_tol` — absolute tolerance.
- `rel_tol` — relative tolerance.
- `why` — explanation of the right answer.

**Outputs**

- `result` — updated test result structure.

## References

- [tests/lib/test_close.m on GitHub (Spinach repository)](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_close.m)
