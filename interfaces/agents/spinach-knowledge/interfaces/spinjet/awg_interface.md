# interfaces/spinjet/awg_interface.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/spinjet/awg_interface.m`
- Signature: `data=awg_interface(spin_system,awg_cmd,cmd_input)`
- Total lines: 141

## Purpose

Interface to the Bruker SpinJet AWG, calling a library of Python scripts that inteface with the Bruker Xepr API. Syntax: data=awg_interface(spin_system,awg_cmd,cmd_input)

## Physical / mathematical content

- This file belongs to the `interfaces` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- spin_system -should contain spinach root directory
- in spin_system.sys.root_dir
- awg_cmd -switch between available commands:
- 'reset_pspel'
- 'compile_pspel_shp'
- 'compile_pspel_def'
- 'compile_pspel_exp'
- 'modify_pspel_defs'
- 'acquire_data'
- cmd_input -cell array of inputs to be passed to
- python scripts

## Outputs

- data.X -abscissa values from the acquired
- signal
- data.rY -real ordinate values from the ac-
- quired signal
- data.iY -imaginary ordinate values from
- the acquired signal.

## Implementation structure

- Interface to the Bruker SpinJet AWG, calling a library of Python
- scripts that inteface with the Bruker Xepr API. Syntax:
- data=awg_interface(spin_system,awg_cmd,cmd_input)
- spin_system -should contain spinach root directory
- in spin_system.sys.root_dir
- awg_cmd -switch between available commands:
- 'reset_pspel'
- 'compile_pspel_shp'
- 'compile_pspel_def'
- 'compile_pspel_exp'
- 'modify_pspel_defs'
- 'acquire_data'

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`, `py_run()`, `report()`, `continue()`, `dlmread()`, `delete()`, `isfield()`, `isfolder()`, `ischar()`, `ismember()`, `iscell()`, `dir()`, `disk()`, `int2str()`.
