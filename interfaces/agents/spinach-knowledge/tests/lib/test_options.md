# tests/lib/test_options.m

- Signature: `options=test_options(varargin)`

## Purpose

Parses name-value options used by the Spinach test runner.

## Parameters / inputs

- `varargin` - name-value pairs. Supported names are `pattern`, `verbose`, and `stop_on_fail`.

## Outputs

- `options` - structure with fields `pattern`, `verbose`, and `stop_on_fail`; defaults are `''`, `false`, and `false`, respectively.

## Implementation structure

- Requires an even number of arguments and rejects unknown option names. `verbose` and `stop_on_fail` accept a logical scalar or the character/string values `true` and `false`; invalid values raise an error. `pattern` is copied from its paired value.
