# tests/lib/new_test_result.m

## Purpose

Creates a regression test result structure for the Spinach test suite. Source: [GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/new_test_result.m).

## Behaviour

- Syntax: `result=new_test_result(id,name,purpose)`.
- Validates the three input arguments via an internal `grumble` helper before building the structure.
- Validation rules:
  - `id` must be a non-empty character row vector, otherwise errors with `'id must be a non-empty character string.'`.
  - `name` must be a non-empty character row vector, otherwise errors with `'name must be a non-empty character string.'`.
  - `purpose` must be a character row vector (empty allowed), otherwise errors with `'purpose must be a character string.'`.
- On success, returns a structure with fields:
  - `id`, `name`, `purpose` — copied from the inputs.
  - `status` — initialised to `'RUNNING'`.
  - `elapsed` — initialised to `0`.
  - `messages` — initialised to `{}`; accumulates one line per check.
  - `failures` — initialised to `{}`; accumulates details of checks that did not pass; an empty `failures` field means no failed check has been recorded.
  - `error` — initialised to `''`.

## Inputs and outputs

Inputs:

- `id` — stable test identifier (non-empty character row vector).
- `name` — short human-readable test name (non-empty character row vector).
- `purpose` — one-sentence purpose statement (character row vector; may be empty).

Output:

- `result` — test result structure as described above.

## References

- Source file: [tests/lib/new_test_result.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/new_test_result.m)
