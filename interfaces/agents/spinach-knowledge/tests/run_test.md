# tests/run_test.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/run_test.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/run_test.m)

## Purpose

Runs a single Spinach regression test selected by an identifier substring, returning the result structure for that test.

## Behaviour

- Validates `test_id` with an internal consistency check (`grumble`), which errors unless the argument is a non-empty character row vector.
- Adds the test library directory (`lib` under the folder containing `run_test.m`) to the MATLAB path.
- Obtains the test manifest via `test_manifest()` and matches `test_id` against both test `id` and `name` fields using substring containment.
- Requires exactly one manifest match; otherwise raises the error `'test_id must match exactly one test.'`.
- Executes the matched test by calling `run_tests` with options `'pattern', test_id, 'verbose', true, 'stop_on_fail', true`.
- Returns the structure produced by `run_tests`; because `run_tests` is invoked with `stop_on_fail` enabled, it throws an error when the matched test fails, so `result` is effectively only returned for a passing test.

## Inputs and outputs

**Syntax**

```
result = run_test(test_id)
```

**Inputs**

- `test_id` — test identifier or unique substring from `list_tests()`; must be a non-empty character row string.

**Outputs**

- `result` — single test result structure; only returned for a passing test, because `run_tests()` throws an error when the matched test fails.

## References

- Spinach GitHub repository: [https://github.com/IlyaKuprov/Spinach](https://github.com/IlyaKuprov/Spinach)
