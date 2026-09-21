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
