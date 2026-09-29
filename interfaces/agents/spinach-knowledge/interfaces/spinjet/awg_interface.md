# interfaces/spinjet/awg_interface.m

[MATLAB implementation](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/spinjet/awg_interface.m) · [Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=awg_interface.m)

## Interface and commands

`data=awg_interface(spin_system,awg_cmd,cmd_input)` dispatches to Python scripts that use the Bruker Xepr API. `spin_system.sys.root_dir` must name an existing directory. `awg_cmd` is a character array and must be one of: `reset_pspel`, `compile_pspel_shp`, `compile_pspel_def`, `compile_pspel_exp`, `modify_pspel_defs`, or `acquire_data`. Omitted `cmd_input` defaults to `{}`; it must be a cell array. The wrapper initialises `data=[]` and only fills it for acquisition.

- `reset_pspel` runs `Xepr_resetexpt`, asks the operator to hide and then show the PulseSPEL window, and pauses for a keypress (Ctrl+C can end the simulation).
- `compile_pspel_shp`, `compile_pspel_def`, and `compile_pspel_exp` pass `cmd_input` to `Xepr_plsspel_shpfile`, `Xepr_plsspel_deffile`, and `Xepr_plsspel_expfile`, respectively. Shape compilation requires exactly one input, the shape-file path, and rejects files at or above 262,144 bytes.
- `modify_pspel_defs` passes the inputs to `Xepr_plsspel_moddefs`.
- `acquire_data` requires exactly three inputs and passes them to `Xepr_getdata`. If the third input is not a character array, it is converted with `num2str`; that value is the temporary-file prefix. The wrapper reads `<prefix>_X.txt`, `<prefix>_rY.txt`, and `<prefix>_iY.txt` with `dlmread` into `data.X`, `data.rY`, and `data.iY`, then deletes those three files.

The source labels `X` as the x-axis and `rY`/`iY` as the real/imaginary signal ordinates, but it does not specify their physical units or impose matrix dimensions; the loaded numeric shapes are those returned by `dlmread`. Other commands leave `data` empty. The script-specific argument meanings belong to the called Python scripts, not this dispatcher.
