# interfaces/spinjet/py_run.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/spinjet/py_run.m`
- Signature: `arg_out=py_run(spin_system,pyscript,arg_in)`
- Total lines: 126

## Purpose

Runs a python script from /interfaces/spinjet/Xepr_python/ folder of Bruker Xepr installation. Variable inputs and output can be passed to and from the script. Syntax: arg_out=py_run(spin_system,pyscript,arg_in)

## Physical / mathematical content

- This file belongs to the `interfaces` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 33-34: Initialise arg_in if not an input; implemented by `if ~exist('arg_in','var'), arg_in={}; end`.
- Lines 36-37: Check consistency; implemented by `grumble(spin_system,pyscript,arg_in)`.
- Lines 39-40: Initialise output as an empty cell; implemented by `arg_out={}`.
- Lines 42-45: Initialise command string to run the selected python script file; implemented by `command_string=['python "' spin_system.sys.root_dir filesep 'interfaces' filesep 'spinjet' filesep 'Xepr_python' filesep pyscript '.py"']`.
- Lines 47-48: Process extra argument cell, if needed; implemented by `if ~isempty(arg_in)`.
- Lines 50-51: Append command string with space; implemented by `command_string=[command_string, ' ']`.
- Lines 53-54: Process each string in the cell of inputs; implemented by `for inputs=1:numel(arg_in)`.
- Lines 56-57: Ensure input as characters; implemented by `if ~ischar(arg_in{inputs})`.
- Lines 61-62: Strip pre-existing surrounding quotes; implemented by `if (numel(arg_in{inputs})>1)&&(arg_in{inputs}(1)=='"')&&(arg_in{inputs}(end)=='"')`.
- Lines 66-67: Append the command string; implemented by `command_string=[command_string '"' arg_in{inputs} '" ']`.
- Lines 73-74: Run Xepr python script with system command-line; implemented by `[status, output_string]=system(command_string)`.
- Lines 76-78: Detect error exception within python script -if no error, then compile outputs from python script, if any were assigned during python print.; implemented by `if ~(status==0)`.

### Control flow inferred from the code

- Line 34: conditional branch on `~exist('arg_in','var'), arg_in={}; end`.
- Line 48: conditional branch on `~isempty(arg_in)`.
- Line 54: `for` loop over `inputs=1:numel(arg_in)`.
- Line 57: conditional branch on `~ischar(arg_in{inputs})`.
- Line 62: conditional branch on `(numel(arg_in{inputs})>1)&&(arg_in{inputs}(1)=='"')&&(arg_in{inputs}(end)=='"')`.
- Line 78: conditional branch on `~(status==0)`.

### Key state/data transformations

- Lines 40: computes `arg_out` using `arg_out={}`.
- Lines 43-45: computes `command_string` using `command_string=['python "' spin_system.sys.root_dir filesep 'interfaces' filesep 'spinjet' filesep 'Xepr_python' filesep pyscript '.py"']`.
- Lines 58: computes `arg_in{inputs}` using `arg_in{inputs}=num2str(arg_in{inputs})`.
- Lines 74: computes `[status, output_string]` using `[status, output_string]=system(command_string)`.

### Local helper functions

- Line 88: `grumble()` — `function grumble(spin_system,py_script_str,python_inputs)`.
  - Representative operation: `if isempty(spin_system)||(~isfield(spin_system,'sys'))|| (~isfield(spin_system.sys,'root_dir'))||(~ischar(spin_system.sys.root_dir))`.
  - Representative operation: `(~isfield(spin_system.sys,'root_dir'))||(~ischar(spin_system.sys.root_dir))`.

## Parameters / inputs

- spin_system -spin_system created by spinach, with
- spin_system.sys.root_dir defined
- pyscript -name of the python script to run, in
- the form of a string without the .py
- extension
- arg_in -cell array containing the inputs to
- the script

## Outputs

- arg_out -cell array of outputs from the script,
- split at every space character (only
- if no error exception has occourred -
- in that case the error message would
- be returned)

## Implementation structure

- Runs a python script from /interfaces/spinjet/Xepr_python/ folder of
- Bruker Xepr installation. Variable inputs and output can be passed to
- and from the script. Syntax:
- arg_out=py_run(spin_system,pyscript,arg_in)
- spin_system -spin_system created by spinach, with
- spin_system.sys.root_dir defined
- pyscript -name of the python script to run, in
- the form of a string without the .py
- extension
- arg_in -cell array containing the inputs to
- the script
- arg_out -cell array of outputs from the script,

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`, `ischar()`, `num2str()`, `system()`, `code()`, `int2str()`, `strsplit()`, `strtrim()`, `isfield()`, `iscell()`.
