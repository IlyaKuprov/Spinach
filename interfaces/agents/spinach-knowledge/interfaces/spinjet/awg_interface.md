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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 39-40: Bootstrap cmd_input if not supplied; implemented by `if ~exist('cmd_input','var'), cmd_input={}; end`.
- Lines 42-43: Check consistency; implemented by `grumble(spin_system,awg_cmd,cmd_input)`.
- Lines 45-46: Initialise data structure; implemented by `data=[]`.
- Lines 48-49: Switch between defined routines; implemented by `switch awg_cmd`.
- Lines 53-54: Reset pulse-spel by connecting with Xepr; implemented by `py_run(spin_system,'Xepr_resetexpt')`.
- Lines 56-57: Instruct user to hide then show the PulseSPEL window (Bruker bug); implemented by `report(spin_system,'<<user: 1) button press hide PulseSPEL window>>')`.
- Lines 64-65: Compile the shape file defined in cmd_input; implemented by `py_run(spin_system,'Xepr_plsspel_shpfile',cmd_input)`.
- Lines 69-70: Compile the definitions file defined in cmd_input; implemented by `py_run(spin_system,'Xepr_plsspel_deffile',cmd_input)`.
- Lines 74-75: Compile the experiment file defined in cmd_input; implemented by `py_run(spin_system,'Xepr_plsspel_expfile',cmd_input)`.
- Lines 79-80: Run command to modify definitions in definitions file; implemented by `py_run(spin_system,'Xepr_plsspel_moddefs',cmd_input)`.
- Lines 84-85: Acquire data from spectometer; implemented by `py_run(spin_system,'Xepr_getdata',cmd_input)`.
- Lines 87-88: Save x axis, imaginary, and real parts of the signal; implemented by `data.X = dlmread([cmd_input{3} '_X.txt'])`.
- Lines 92-93: Delete temporary file of results; implemented by `delete([cmd_input{3} '_X.txt'],[cmd_input{3} '_rY.txt'],[cmd_input{3} '_iY.txt'])`.

### Control flow inferred from the code

- Line 40: conditional branch on `~exist('cmd_input','var'), cmd_input={}; end`.
- Line 49: dispatches on `awg_cmd`; cases `{'reset_pspel'}`, `{'compile_pspel_shp'}`, `{'compile_pspel_def'}`, `{'compile_pspel_exp'}`, `{'modify_pspel_defs'}`, `{'acquire_data'}`.

### Key state/data transformations

- Lines 46: computes `data` using `data=[]`.
- Lines 88: computes `data.X` using `data.X = dlmread([cmd_input{3} '_X.txt'])`.
- Lines 89: computes `data.rY` using `data.rY = dlmread([cmd_input{3} '_rY.txt'])`.
- Lines 90: computes `data.iY` using `data.iY = dlmread([cmd_input{3} '_iY.txt'])`.

### Local helper functions

- Line 100: `grumble()` — `function grumble(spin_system,awg_cmd,cmd_input)`.
  - Representative operation: `if ~isfield(spin_system,'sys') || ~isfield(spin_system.sys,'root_dir')`.
  - Representative operation: `error('spin_system.sys.root_dir must be supplied')`.

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
