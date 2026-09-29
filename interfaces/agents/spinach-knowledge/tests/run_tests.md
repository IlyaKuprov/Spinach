# tests/run_tests.m

## Purpose

Runs the Spinach regression test suite ([source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/run_tests.m)).

## Behaviour

- Syntax: `results=run_tests(varargin)`.
- Enforces that options are supplied as name-value pairs; option names must be non-empty character row strings, otherwise an error is raised.
- Adds the test library to the path: the directory of the file itself, its `lib` subdirectory, and recursively `kernel` and `interfaces` subdirectories.
- Adds the Spinach production directories to the path: recursively `etc`, `experiments`, `interfaces`, and `kernel` under the Spinach root (the parent of the test directory).
- Parses options via `test_options`, obtains the test manifest via `test_manifest`.
- If a `pattern` option is non-empty, filters the manifest to tests whose `id` or `name` contains the pattern as a substring.
- Preallocates an empty result structure array with fields `id`, `name`, `purpose`, `status`, `elapsed`, `messages`, `failures`, `error`.
- For each manifest entry: resets a clock and the record accumulator via `test_record(new_test_result(id,name,''))`, then runs the test function with `feval`.
- On normal return, sets `elapsed` from `toc`, sets `error` to `strjoin(result.failures,'; ')`, and sets `status` to `'PASS'` if `failures` is empty, otherwise `'FAIL'`.
- On a MATLAB error, recovers the recorded result via `test_record([])`, sets `status` to `'FAIL'`, sets `elapsed` from `toc`, appends `err.message` to `failures`, and sets `error` to `strjoin(result.failures,'; ')`.
- A test is reported as a failure when it records at least one failed check in the `failures` field of its result structure, in which case the `error` field carries the recorded failures; messages accumulated before the failure are preserved. A test that throws a MATLAB error is also reported as a failure, and the suite continues with the next test; checks recorded before the error are recovered from `test_record` and reported alongside the error message.
- If `verbose` is set, prints per-test status, id, and elapsed time in the format `%s\t%s\t%.3f s`, followed by indented messages, and, for failures, an `ERROR:` line with the joined error string.
- If `stop_on_fail` is set and a test fails, the loop breaks and no further tests run.
- After the loop, counts `'PASS'` and `'FAIL'` statuses and prints `Spinach regression tests: %d passed, %d failed.`.
- If any test failed, prints a `FAIL\t%s\t%s` line per failed test with its id and error string, then raises the error `'Spinach regression tests failed.'`.

## Inputs and outputs

Inputs:

- `varargin` - name-value options: `'pattern'`, `'verbose'`, and `'stop_on_fail'`.

Outputs:

- `results` - structure array with test outcomes and messages, with fields `id`, `name`, `purpose`, `status`, `elapsed`, `messages`, `failures`, `error`.

## References

- [Source code](https://github.com/IlyaKuprov/Spinach/blob/main/tests/run_tests.m)
