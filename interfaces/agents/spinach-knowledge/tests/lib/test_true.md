# tests/lib/test_true.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_true.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_true.m)

## Purpose

`test_true` adds a logical regression check with a clear message to a Spinach test suite. It evaluates a pass/fail condition, appends a human-readable PASS or FAIL message to the running test result structure, and returns a `passed` flag so the caller can guard subsequent statements.

## Behaviour

- Syntax: `[result,passed]=test_true(result,label,condition,why)`.
- Input consistency is enforced first by the local `grumble` function, which errors when:
  - `result` is not a scalar structure with `messages` and `failures` fields (`'result must be a scalar test result structure.'`),
  - `label` is not a non-empty row character string (`'label must be a non-empty character string.'`),
  - `condition` is neither logical nor numeric (`'condition must be logical or numeric.'`),
  - `why` is not a non-empty row character string (`'why must be a non-empty character string.'`).
- The condition is evaluated as `passed=isscalar(condition)&&(condition~=0)`; a check passes only when `condition` is scalar and non-zero.
- On pass, the string `'PASS: ' label ' -- ' why` is appended to `result.messages`.
- On failure, `'FAIL: ' label ' -- ' why` is appended to `result.messages` and `label ' -- ' why` is also appended to `result.failures`.
- A failed check is recorded in the `messages` and `failures` fields of the result structure rather than thrown, so that the checks that follow it in the calling test are still evaluated and the earlier passes are not lost; `run_tests` inspects the `failures` field to decide the status.
- Where the check is a precondition for the statements that follow it, the caller must guard those statements with the `passed` flag, because a secondary Matlab error would send the test into the catch path of `run_tests` and leave the remaining checks unevaluated.
- After every check, `test_record(result)` is called so the record is retained; the catch path of `run_tests` still reports the checks completed before such an error.

## Inputs and outputs

**Inputs**

- `result` — test result structure.
- `label` — check label.
- `condition` — logical pass/fail condition (logical or numeric).
- `why` — explanation of the right answer.

**Outputs**

- `result` — updated test result structure.
- `passed` — true when the condition held.

## References

- Source file: [tests/lib/test_true.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_true.m) in the Spinach repository on GitHub.
- Related functions referenced in the source comments: `run_tests`, `test_record`.
