# tests/lib/test_close.m

- Signature: `result=test_close(result,label,observed,reference,abs_tol,rel_tol,why)`

## Purpose

Compares an observed array with a reference array using absolute and relative tolerances, and records the check in the test result.

## Numerical / algorithmic content

- Sparse inputs are converted to full arrays and the values to double vectors. Arrays must have matching dimensions and finite values.
- The error is the Euclidean norm of the vector difference. The allowed limit is `abs_tol+rel_tol*max(1,norm(reference,2))`. A size mismatch, non-finite input or comparison, or error above the limit is recorded as a failure.

## Parameters / inputs

- `result` - scalar test-result structure with `messages` and `failures` fields.
- `label` - non-empty character-row check label.
- `observed`, `reference` - numeric or logical arrays to compare.
- `abs_tol`, `rel_tol` - non-negative real numeric scalar tolerances.
- `why` - non-empty character-row explanation included in the result message.

## Outputs

- `result` - updated with a PASS or FAIL message; failures are appended to `failures`.

## Implementation structure

- Validates inputs, records the comparison outcome, then calls `test_record` so the completed check is retained for the test runner.
