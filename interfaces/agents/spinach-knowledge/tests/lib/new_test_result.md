# tests/lib/new_test_result.m

- Signature: `result=new_test_result(id,name,purpose)`

## Purpose

Initialises a regression-test result structure for subsequent checks.

## Parameters / inputs

- `id` - non-empty character-row test identifier.
- `name` - non-empty character-row test name.
- `purpose` - character-row purpose text; it may be empty.

## Outputs

- `result` - structure containing the supplied `id`, `name`, and `purpose`; status `RUNNING`; elapsed time `0`; empty `messages` and `failures` cells; and an empty `error` string.

## Implementation structure

- Validates the three inputs, then initializes the result fields. Invalid `id` or `name` values, or a non-character or non-row `purpose`, raise an error.
