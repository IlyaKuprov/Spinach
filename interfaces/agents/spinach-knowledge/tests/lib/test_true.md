# tests/lib/test_true.m

- Signature: `result=test_true(result,label,condition,why)`

## Purpose

Records a logical regression check and returns its updated test-result structure. The function also provides a `passed` output when requested.

## Parameters / inputs

- `result` - scalar test-result structure with `messages` and `failures` fields.
- `label` - non-empty character-row check label.
- `condition` - logical or numeric condition; only a nonzero scalar passes.
- `why` - non-empty character-row explanation included in the result message.

## Outputs

- `result` - updated with a PASS or FAIL message; failures are appended to `failures`.
- `passed` - optional logical scalar indicating whether the condition passed.

## Implementation structure

- Validates the inputs, evaluates whether `condition` is a nonzero scalar, records the outcome, and calls `test_record` for the runner's catch path.
