# fix_path.m

- Signature: `fix_path(config_style)`

## Purpose

Configure MATLAB's search path for Spinach. With no argument, the function uses the `noob` configuration. It can also add Spinach to the current path, remove Spinach folders, or reset MATLAB's path before adding Spinach.

## Physical / mathematical content

Not applicable: this utility changes the MATLAB search path and does not perform a physical or mathematical calculation.

## Numerical / algorithmic content

No numerical algorithm is used. The selected configuration determines which path operations are performed.

## Implementation structure

- Default `config_style` to `noob` when the argument is omitted.
- Check that `config_style` is a character string; reject other types.
- Locate the Spinach root from this file's full path.
- For `noob` or `reset`, restore MATLAB's default path, add the `etc`, `experiments`, `interfaces`, and `kernel` trees, then run `existentials`.
- For `add`, add those four trees to the existing path and run `existentials`.
- For `remove`, remove those four trees from the path.
- Report the operation in the MATLAB console; reject any unrecognized configuration with an error.
