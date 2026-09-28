# interfaces/spinjet/awg_interface.m

- Signature: `data=awg_interface(spin_system,awg_cmd,cmd_input)`

## Purpose

Interface to the Bruker SpinJet AWG through Python scripts that interact with the Bruker Xepr API.

## Parameters / inputs

- `spin_system` must contain `spin_system.sys.root_dir`, a path to an existing directory.
- `awg_cmd` must be one of the following character-array commands:
  - `'reset_pspel'`: reset PulseSPEL through Xepr.
  - `'compile_pspel_shp'`: compile the shape file specified by `cmd_input`.
  - `'compile_pspel_def'`: compile the definitions file specified by `cmd_input`.
  - `'compile_pspel_exp'`: compile the experiment file specified by `cmd_input`.
  - `'modify_pspel_defs'`: modify definitions in the definitions file.
  - `'acquire_data'`: acquire data from the spectrometer.
- `cmd_input` is a cell array of inputs passed to the Python scripts. If omitted, it defaults to `{}`. The three compile commands each require exactly one element. For `'compile_pspel_shp'`, the specified shape file must be smaller than 262144 bytes on disk. `'modify_pspel_defs'` requires an even number of elements, at least two. `'acquire_data'` requires exactly three elements; its third element is used as the result-file prefix and is converted to a string with `num2str` if it is not a character array.

## Outputs

- `data.X`: abscissa values of the acquired signal.
- `data.rY`: real ordinate values of the acquired signal.
- `data.iY`: imaginary ordinate values of the acquired signal.

`data` is initialized to `[]`; these fields are populated for `'acquire_data'`.

## Implementation structure

Each command invokes its corresponding Xepr Python script. After `'reset_pspel'`, the function instructs the user to hide and then show the PulseSPEL window and pauses for a keypress. For `'acquire_data'`, it reads `<prefix>_X.txt`, `<prefix>_rY.txt`, and `<prefix>_iY.txt` into the corresponding output fields, then deletes those temporary files.