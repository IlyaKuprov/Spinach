# tests/lib/test_options.m

## Purpose

Parses name-value options for the Spinach test runner.

## Behaviour

- Sets defaults: `pattern` is an empty character array, `verbose` is `false`, and `stop_on_fail` is `false`.
- Errors with `'options must be supplied as name-value pairs.'` if the number of arguments is not even.
- Iterates over the arguments in steps of two, treating each even-indexed argument as an option name and the following argument as its value.
- For `'pattern'`, assigns the following value directly to `options.pattern`.
- For `'verbose'` and `'stop_on_fail'`, accepts a scalar logical, or the character or scalar string `'true'`/`'false'`; any other value raises an error (`'verbose option must be true or false.'` or `'stop_on_fail option must be true or false.'`).
- Any unrecognised option name raises an error `'unknown option: ' varargin{n}`.

## Inputs and outputs

- `varargin` — name-value option pairs.
- `options` — options structure with fields `pattern`, `verbose`, and `stop_on_fail`.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_options.m)
