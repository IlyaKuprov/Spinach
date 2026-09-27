# interfaces/spinjet/py_run.m

- Signature: `arg_out=py_run(spin_system,pyscript,arg_in)`

## Purpose

Runs a Python script from the `interfaces/spinjet/Xepr_python/` folder of a Bruker Xepr installation, passing inputs to the script and returning its printed output.

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- `spin_system` - spin system created by Spinach, with `spin_system.sys.root_dir` defined as a character string.
- `pyscript` - name of the Python script as a character string, without the `.py` extension.
- `arg_in` - cell array containing inputs to the script. If omitted, it defaults to an empty cell array. Non-character entries are converted with `num2str`.

## Outputs

- `arg_out` - cell array of tokens from the script's standard output, split on whitespace. It is empty if the command succeeds without producing output. If the command returns a nonzero status, the function raises an error containing the status code and returned output rather than returning `arg_out`.

## Implementation structure

The function checks the spin system, script name, script file, and input cell array before running the script. It builds a `python` command using `spin_system.sys.root_dir` and the script name. Each supplied input is passed as a double-quoted command-line argument; existing surrounding double quotes are removed, and backslashes, dollar signs, backticks, and double quotes are escaped. It then runs the command with `system` and splits nonempty output after trimming it.