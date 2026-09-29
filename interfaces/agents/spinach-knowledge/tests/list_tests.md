# tests/list_tests.m

## Purpose

Lists the Spinach regression tests available in the `tests` directory ([source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/list_tests.m)).

## Behaviour

The function adds the test library (`lib` subdirectory of the folder containing `list_tests.m`) to the MATLAB path, parses the supplied options, and obtains the test manifest from `test_manifest()`. If a non-empty `pattern` option is given, the manifest is filtered to the entries whose `id` or `name` contains the pattern as a substring (case-sensitive `contains` check). The remaining entries are printed to the command window, one per line, as identifier, a tab, and the test name.

## Inputs and outputs

**Syntax**

```matlab
manifest = list_tests(varargin)
```

**Inputs**

- `varargin` — optional name-value pair `'pattern'`, a string used as a substring filter on test identifiers and names; parsed by `test_options`.

**Outputs**

- `manifest` — structure array with test identifiers (`id`) and names (`name`), as returned by `test_manifest()` and optionally filtered by the pattern.

## References

- [Spinach regression test lister — source file](https://github.com/IlyaKuprov/Spinach/blob/main/tests/list_tests.m)
