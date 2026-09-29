# interfaces/spinjet/py_run.m

[MATLAB implementation](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/spinjet/py_run.m) · [Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=py_run.m)

## Interface

`arg_out=py_run(spin_system,pyscript,arg_in)` runs `<spin_system.sys.root_dir>/interfaces/spinjet/Xepr_python/<pyscript>.py` using the `python` command and MATLAB `system`. `spin_system.sys.root_dir` must be present as a character array, and `pyscript` must be a non-empty character array; the named script file must exist. The `.py` suffix is added by the wrapper. If `arg_in` is omitted it defaults to `{}`; supplied inputs must be a cell array.

Each cell input is converted with `num2str` unless it is already a character array. The wrapper strips one pair of surrounding double quotes when present, escapes backslashes, dollar signs, backticks, and double quotes, and appends each value as a double-quoted command-line argument. It does not marshal MATLAB arrays as Python objects: the script receives textual command-line arguments.

A nonzero `system` status raises an error containing the status and Python output. On success, the wrapper trims the captured stdout and splits it into whitespace-delimited text tokens, returning those tokens as a cell array (transposed from `strsplit`'s row result). It does not parse the tokens into numeric or structured MATLAB values. No cache or physical-unit semantics are defined in this helper.
