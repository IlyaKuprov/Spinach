# kernel/overloads/save_anyway.m

Source: [kernel/overloads/save_anyway.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/save_anyway.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=save_anyway.m)

- Signature: `save_anyway(file_name,variable)`

## Purpose

A wrapper intended to let SPMD blocks save a value using a fixed variable name in a MAT-file.

## Behaviour

The function checks only that `file_name` is a character array, then calls `save(file_name,'variable','-v7.3')`. Consequently, the saved MAT-file variable is named `variable`, regardless of the caller's name for the second argument. It saves one variable per call in MATLAB v7.3 format and then calls `drawnow`. There is no returned output; save errors are not caught by this wrapper.

## Inputs

- `file_name`: character array naming the MAT-file.
- `variable`: value to save.

The source does not add validation for the value being saved.
