# kernel/utilities/idxof.m

## Purpose

`idxof.m` returns the index of a spin in the Spinach input structure given its label, allowing interactions to be specified by spin label rather than by number.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/idxof.m>

## Behaviour

- Syntax: `idx=idxof(sys,label)`.
- The function first runs a consistency check (`grumble`) on the inputs.
- It locates the label by scanning `sys.labels` with `strcmp` via `cellfun` and `find`.
- If no match is found, it errors with `'label not found.'`.
- If more than one match is found, it errors with `'labels are not unique'`.
- The consistency check errors when: `label` is not a character array (`'label must be a character string.'`), `sys.labels` is missing (`'sys.labels is missing.'`), or the number of elements in `sys.labels` does not match the number of elements in `sys.isotopes` (`'number of elements in sys.label must match the number of spins.'`).

## Inputs and outputs

Inputs:

- `sys` — Spinach input structure that includes a `sys.labels` field with unique labels.
- `label` — label whose index is to be returned; must be a character string.

Outputs:

- `idx` — the index of the spin, an integer.

## References

- Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=idxof.m>
- Source file: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/idxof.m>
