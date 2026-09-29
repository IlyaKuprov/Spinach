# fix_path.m

- Signature: `fix_path(config_style)`
- Return value: none; the function changes MATLAB's search path and reports status to the command window.

## Purpose

Sets up or removes Spinach directories on the active MATLAB path. The Spinach root is located relative to `fix_path.m`; the four managed trees are `etc`, `experiments`, `interfaces`, and `kernel`, including their subdirectories.

## Accepted input

`config_style` must be a character array. If omitted, it defaults to `'noob'`. The accepted values and effects are:

- `'noob'` and `'reset'`: call MATLAB's `restoredefaultpath`, add the four Spinach trees to the beginning of the path, then run `existentials` checks. This resets the MATLAB path before adding Spinach.
- `'add'`: preserve the existing MATLAB path, add those same Spinach trees at the beginning, then run `existentials` checks.
- `'remove'`: remove those Spinach trees from the path and report the removal; it does not reset MATLAB's path or remove unrelated entries.

An unrecognised style raises an error. Because the type check uses `ischar`, a MATLAB string scalar is not the documented character-array input.

## Operational notes

The path edits use `genpath` over each managed tree. The reset modes therefore replace the current path with MATLAB's default path before adding Spinach; use `'add'` when unrelated existing path entries should be retained. No value is returned.

## Links

- MATLAB source: https://github.com/IlyaKuprov/Spinach/blob/main/fix_path.m
- [Spinach Wiki: fix_path.m](https://spindynamics.org/wiki/index.php?title=fix_path.m)
