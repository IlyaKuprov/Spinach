# kernel/utilities/report.m

- Signature: `report(spin_system,report_string)`

## Purpose

Prints a character message with a prefix derived from the call stack. A one-argument call errors; an empty `spin_system` selects standard output. Setting `spin_system.sys.output='hush'` suppresses the call, and otherwise the message is written with a trailing newline to the configured output destination.

## Physical / mathematical content

This is an output utility; it does not perform spin dynamics or numerical calculations.

## Numerical / algorithmic content

The function removes known uninformative stack entries, labels parallel-worker stack entries as `parfor/spmd`, reverses and joins the remaining caller names, removes the final three-character suffix, and pads the prefix to 50 characters. Longer prefixes are shortened to 50 characters with a leading ellipsis. It then prepends the prefix and writes the result with `fprintf`; impossible write errors are suppressed.

## Parameters / inputs

- `spin_system` - system structure with `sys.output` set to `'hush'` or a file identifier; an empty value defaults to output identifier 1.
- `report_string` - character array to report; non-character input errors.

## Outputs

No value is returned. The message is printed to the configured destination, unless output is hushed.

## Implementation structure

Validation requires `spin_system.sys.output` and a character `report_string`. With output enabled, the function builds the call-stack prefix and writes `[prefix ]  report_string` followed by a newline. Errors from the final `fprintf` are caught and ignored.

## Reference

[Spin Dynamics Wiki: report.m](https://spindynamics.org/wiki/index.php?title=report.m)
