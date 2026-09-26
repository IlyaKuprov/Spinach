# interfaces/spinjet/awg_interface.m

- Signature: `data=awg_interface(spin_system,awg_cmd,cmd_input)`

## Purpose

Interface to the Bruker SpinJet AWG, calling a library of Python scripts that inteface with the Bruker Xepr API. Syntax: data=awg_interface(spin_system,awg_cmd,cmd_input)

## Physical / mathematical content

## Numerical / algorithmic content

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
