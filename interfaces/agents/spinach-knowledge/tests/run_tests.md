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
- Each test runs through `test_execute`, which captures console output while preserving recorded failures, partial check messages, and thrown-error text. Elapsed time uses a dedicated timer so a test's own `tic` calls cannot reset it.
- Progress includes the current index, selected count, test identifier, status, and elapsed time. `verbose=true` additionally prints recorded check explanations. Captured chatter is printed only for failures, including output preceding exceptions.
- Test graphics are invisible. Cleanup deletes newly created figures (including hidden handles), restores root graphics defaults and the previously current figure, and retains existing figures. Tests must not explicitly delete caller-owned figures.
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
