# interfaces/spinjet/py_run.m

- Signature: `arg_out=py_run(spin_system,pyscript,arg_in)`

## Purpose

Runs a python script from /interfaces/spinjet/Xepr_python/ folder of Bruker Xepr installation. Variable inputs and output can be passed to and from the script. Syntax: arg_out=py_run(spin_system,pyscript,arg_in)

## Physical / mathematical content

## Numerical / algorithmic content

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
