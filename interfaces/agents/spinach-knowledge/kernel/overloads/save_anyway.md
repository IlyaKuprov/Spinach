# kernel/overloads/save_anyway.m

- Signature: `save_anyway(file_name,variable)`

## Purpose

A wrapper intended to trick SPMD blocks into saving data. Can only save one variable at a time, its name in the mat file is "variable". Syntax: save_anyway(file_name,variable)

## Physical / mathematical content

None; this is a file-saving utility.

## Numerical / algorithmic content

The function checks that `file_name` is a character array, saves the input as the MAT-file variable `variable` using `-v7.3`, and calls `drawnow`.

## Parameters / inputs

- `file_name` — character string specifying the file name.
- `variable` — the variable to save.

## Outputs

No output arguments.
