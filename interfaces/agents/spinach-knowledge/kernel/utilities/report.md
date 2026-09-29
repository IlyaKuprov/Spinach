# kernel/utilities/report.m

## Purpose

Writes a log message to the console or an ASCII file, prefixed with the call stack of the function that produced it. Source: [kernel/utilities/report.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/report.m).

## Behaviour

- Syntax: `report(spin_system,report_string)`.
- A single-argument call raises the error `console reporting function requires two arguments.`.
- If `spin_system` is empty, `spin_system.sys.output` is set to `1` (console) before proceeding.
- If `spin_system.sys.output` is `'hush'`, the call is ignored and nothing is printed.
- Otherwise the input is validated by the internal `grumble` function, which errors when:
  - `spin_system.sys` or `spin_system.sys.output` does not exist: `spin_system.sys.output field must exist.`;
  - `spin_system.sys.output` is neither `double` nor `char`, or is `char` but not exactly `'hush'`: `spin_system.sys.output must be either 'hush' or a file ID.`;
  - `report_string` is not a character array: `report_string must be a string.`.
- The call stack is obtained with `dbstack`. Entries listed as uninformative (`hamiltonian/parfor_progr`, `@(~)parfor_progr`, `iDispatchDataReceived`, `@(src,event)iDispatchDataReceived(func,src,event)`, `DataQueue.dispatchContinuation`, `DataQueue.maybeDrainAndDispatchAllDataOnQueue`, `AbstractDataQueue.notifyQueue`, `ParforEngine.getCompleteIntervals`, `parallel_function`, `powder/parfor_progr`, `make_general_channel/channel_general`, `ThreadsParforEngine.getCompleteIntervals`) are deleted; `distributed_execution` is replaced with `parfor/spmd > `; every other entry is suffixed with ` > `.
- The stack names are concatenated in reverse order (from caller down to `report`), the trailing `.m`-style extension is trimmed by dropping the last three characters, and an empty prefix is replaced by a single space.
- The prefix is rolled to a fixed width: if shorter than 50 characters it is padded to 50; otherwise it is truncated to `'...'` plus the last 47 characters.
- The final line is `'[' prefix ' ]  ' report_string`, written with `fprintf(spin_system.sys.output,'%s\n',...)` inside a `try`/`end` block that ignores impossible writes; a trailing newline is appended automatically, so the input string need not end with one.
- All output can be silenced by setting `sys.output='hush'` in the Spinach input stream or by setting `spin_system.sys.output='hush'` at any point during the calculation.

## Inputs and outputs

**Inputs**

- `spin_system` — Spinach system object; its `sys.output` field must be `'hush'` or a file ID (`double` or `char`). An empty `spin_system` is tolerated and treated as console output.
- `report_string` — character string with the message to log.

**Outputs**

- None returned; the function prints the message to the console or to the destination specified in `spin_system.sys.output`.

## References

1. Spinach Wiki: [report.m](https://spindynamics.org/wiki/index.php?title=report.m)
2. Source file: [kernel/utilities/report.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/report.m)
