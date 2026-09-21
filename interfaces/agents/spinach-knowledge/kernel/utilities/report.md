# kernel/utilities/report.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/report.m`
- Signature: `report(spin_system,report_string)`
- Total lines: 126

## Purpose

Writes a log message to the console or an ACSII file. The message includes the call stack of the function that produced it. Syntax: report(spin_system,report_string)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- report_string -a character string

## Outputs

- this function prints the message to the console or to the
- destination specified in spin_system.sys.output
- Note: a newline symbol at the end of the string is not neces-
- sary -it is added by the function.
- Note: all output produced by this function may be silenced
- by setting sys.output='hush' in the Spinach input
- stream or by setting spin_system.sys.output='hush'
- at any point during the calculation.

## Implementation structure

- Writes a log message to the console or an ACSII file. The message
- includes the call stack of the function that produced it. Syntax:
- report(spin_system,report_string)
- report_string -a character string
- this function prints the message to the console or to the
- destination specified in spin_system.sys.output
- Note: a newline symbol at the end of the string is not neces-
- sary -it is added by the function.
- Note: all output produced by this function may be silenced
- by setting sys.output='hush' in the Spinach input
- stream or by setting spin_system.sys.output='hush'
- at any point during the calculation.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `strcmp()`, `grumble()`, `iDispatchDataReceived()`, `ismember()`, `call_stack()`, `prefix_string()`, `pad()`, `isfield()`, `ischar()`.
