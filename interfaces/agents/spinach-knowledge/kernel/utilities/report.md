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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 30-31: Catch single-argument calls; implemented by `if nargin==1`.
- Lines 35-36: Catch outside calls with empty spin_system object; implemented by `if isempty(spin_system), spin_system.sys.output=1; end`.
- Lines 38-39: Ignore the call if the system is hushed; implemented by `if ~strcmp(spin_system.sys.output,'hush')`.
- Lines 41-42: Validate the input; implemented by `grumble(spin_system,report_string)`.
- Lines 44-50: List call stack exceptions; implemented by `not_useful={'hamiltonian/parfor_progr','@(~)parfor_progr','iDispatchDataReceived', 'iDispatchDataReceived','@(src,event)iDispatchDataReceived(func,src,event)', 'DataQueu…`.
- Lines 53-54: Get the call stack; implemented by `call_stack=dbstack`.
- Lines 56-57: Compose the prefix; implemented by `for n=1:numel(call_stack)`.
- Lines 59-60: Uninformative entries; implemented by `if ismember(call_stack(n).name,not_useful)`.
- Lines 62-63: Simply delete; implemented by `call_stack(n).name=''`.
- Lines 65-66: Parallelisation rigging entries; implemented by `elseif ismember(call_stack(n).name,mdcs)`.
- Lines 68-69: Tell the user it's a parallel call; implemented by `call_stack(n).name='parfor/spmd > '`.
- Lines 73-74: Just add the entry; implemented by `call_stack(n).name=[call_stack(n).name ' > ']`.
- Lines 80-81: Concatenate the names; implemented by `prefix_string=[call_stack(end:-1:2).name]`.
- Lines 83-84: Remove m-file extensions; implemented by `prefix_string=prefix_string(1:(end-3))`.
- Lines 86-87: Fix occasional empty prefixes; implemented by `if isempty(prefix_string), prefix_string=' '; end`.
- Lines 89-90: Roll the prefix; implemented by `if numel(prefix_string)<50`.
- Lines 96-97: Add prefix to the report string; implemented by `report_string=['[' prefix_string ' ] ' report_string]`.
- Lines 99-100: Send the report string to the output, ignoring impossible writes; implemented by `try fprintf(spin_system.sys.output,'%s\n',report_string); end`.

### Control flow inferred from the code

- Line 31: conditional branch on `nargin==1`.
- Line 36: conditional branch on `isempty(spin_system), spin_system.sys.output=1; end`.
- Line 39: conditional branch on `~strcmp(spin_system.sys.output,'hush')`.
- Line 57: `for` loop over `n=1:numel(call_stack)`.
- Line 60: conditional branch on `ismember(call_stack(n).name,not_useful)`.
- Line 87: conditional branch on `isempty(prefix_string), prefix_string=' '; end`.
- Line 90: conditional branch on `numel(prefix_string)<50`.

### Key state/data transformations

- Lines 45-50: computes `not_useful` using `not_useful={'hamiltonian/parfor_progr','@(~)parfor_progr','iDispatchDataReceived', 'iDispatchDataReceived','@(src,event)iDispatchDataReceived(func,src,event)', 'DataQueu…`.
- Lines 51: computes `mdcs` using `mdcs={'distributed_execution'}`.
- Lines 54: computes `call_stack` using `call_stack=dbstack`.
- Lines 63: computes `call_stack(n).name` using `call_stack(n).name=''`.
- Lines 81: computes `prefix_string` using `prefix_string=[call_stack(end:-1:2).name]`.
- Lines 97: computes `report_string` using `report_string=['[' prefix_string ' ] ' report_string]`.

### Local helper functions

- Line 107: `grumble()` — `function grumble(spin_system,report_string)`.
  - Representative operation: `if (~isfield(spin_system,'sys'))||~isfield(spin_system.sys,'output')`.
  - Representative operation: `error('spin_system.sys.output field must exist.')`.

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
